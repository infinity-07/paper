% =========================================================================
% 网格与索引示意图 (周期边界扩展 + WENO 重构)
% 
% 索引说明：
%   |     cell(i-2)     |      cell(i-1)     |       cell(i)       |      cell(i+1)     |      cell(i+2)     |      cell(i+3)     |
%   |                   |                    |                     |                    |                    |                    |
% ---o------------------o--------------------o---------------------o--------------------o--------------------o--------------------o---
%    ^                  ^                    ^                     ^                    ^                    ^                    ^
% boundary(i-5/2)   boundary(i-3/2)     boundary(i-1/2)       boundary(i+1/2)     boundary(i+3/2)     boundary(i+5/2)     boundary(i+7/2)
%
% 在每个boundary上，我们需要重构出 boundary^- 和 boundary^+。
% 
% 在有限体积方法中，数值通量的计算有两种常见的索引方式：
% 
% 1) 以 cell 
%    - cell 的数量为 n
%    - 每个 cell 有左右两个边界, 
%    - 对于 cell(i) 应当重构出 boundary(i-1/2)^+ 和 boundary(i+1/2)^-
%
%
% 2) 以 boundary 为基准：
%    - boundary 的数量为 n+1
%    - 每个 boundary 有左右两个重构值
%    - 对于 boundary(i+1/2) 应当重构出 boundary(i+1/2)^- 和 boundary(i+1/2)^+
%
% 无论采用哪种方式，都需要在边界之外补充额外的点，以满足 WENO 重构 stencil 的需求。
% 
% 个人认为以 cell 为索引更直观，因此本程序采用 cell 索引 的方式。
%
%--------------------------------------------------------------------------
% 下面以 boundary(i+1/2) 为例介绍重构模板
%   - boundary(i+1/2)^- 使用 cell(i-2), cell(i-1), cell(i), cell(i+1), cell(i+2) 重构
%   - boundary(i+1/2)^+ 使用 cell(i-1), cell(i), cell(i+1), cell(i+2), cell(i+3) 重构
% 
% 在本代码中，weno 函数的输入变量为 u(i-2), u(i-1), u(i), u(i+1), u(i+2)，输出变量为 boundary(i+1/2)^-。
% 由于 weno 函数具有一定的对称性，因此对于 boundary(i+1/2)^+，我们只需将输入变量顺序反转即可。(HWENO 不满足该性质)
% 也就是  boundary(i+1/2)^- = weno(u(i-2), u(i-1), u(i), u(i+1), u(i+2))
%         boundary(i+1/2)^+ = weno(u(i+3), u(i+2), u(i+1), u(i), u(i-1))
%
% -------------------------------------------------------------------------
% 本程序实现思路：
%
% - 输入解向量 u (长度为 num_points)
% - 在两端各扩展 3 个 ghost cell，形成 u_extended (长度为 num_points+6)
% - 使用 WENO 重构，在每个 cell 左右边界计算重构值
% - 将结果存入 data 结构体
%
% data 的定义：
%   data(i).leftBoundary   = 第 i 个 cell 左边界的重构值
%   data(i).rightBoundary  = 第 i 个 cell 右边界的重构值
%
% 实际后续计算中只会用到：
%   - data(4 : num_points+3) 的左右重构值
%   - data(3) 的右重构值，以及 data(num_points+4) 的左重构值
% 为了保持索引一致性，这里额外计算了 data(3).leftBoundary 和 data(num_points+4).rightBoundary。
%
% -------------------------------------------------------------------------
% 原始区间:   u(1) ----------- u(num_points)
%
% 周期边界扩展后:
%   u_extended = [ u(end-2)  u(end-1)  u(end) | u(1)  u(2)  u(3) ... u(num_points) | u(1)  u(2)  u(3) ]
%
% =========================================================================

function rhs = getL(u, flux, flux_d, dx)
global config;
global lambdaVec
lambdaVec = zeros(1, 206);

    % 计算数值通量的空间导数 (采用 WENO + LLF)

    % 网格点数
    N = length(u);
    
    % 周期边界扩展：左右各加 3 个 ghost cells
    u_ext = [u(end-2:end), u, u(1:3)];

    % 最大特征速度 (LLF 通量)
    alpha = max(abs(flux_d(u)));
    
    % 初始化边界值
    faces(1:N+6) = struct('left', 0.0, 'right', 0.0);

    if config ~= 4
        % WENO 重构
        for i = 3:N+4
            % 右边界重构值 boundary(i+1/2)^-
            [ur, lambda] = wenoSwitch(u_ext(i-2), u_ext(i-1), u_ext(i), u_ext(i+1), u_ext(i+2));
            lambdaVec(i) = lambdaVec(i)+lambda;

            % 左边界重构值 boundary(i-1/2)^+
            [ul, lambda] = wenoSwitch(u_ext(i+2), u_ext(i+1), u_ext(i), u_ext(i-1), u_ext(i-2));
            lambdaVec(i) = lambdaVec(i)+lambda;

            faces(i).left  = ul;
            faces(i).right = ur;
        end
    else
        facesPol(1:N+6) = struct('left', 0.0, 'right', 0.0);
        facesExp(1:N+6) = struct('left', 0.0, 'right', 0.0);

        for i = 3:N+4
            % 右边界重构值 boundary(i+1/2)^-
            ur = wenoPol3(u_ext(i-1), u_ext(i), u_ext(i+1));
            % 左边界重构值 boundary(i-1/2)^+
            ul = wenoPol3(u_ext(i+1), u_ext(i), u_ext(i-1));

            facesPol(i).left  = ul;
            facesPol(i).right = ur;

            faces(i).left  = ul;
            faces(i).right = ur;
        end

        for i = 3:N+4
            % 右边界重构值 boundary(i+1/2)^-
            ur = wenoExp3(u_ext(i-1), u_ext(i), u_ext(i+1));
            % 左边界重构值 boundary(i-1/2)^+
            ul = wenoExp3(u_ext(i+1), u_ext(i), u_ext(i-1));

            facesExp(i).left  = ul;
            facesExp(i).right = ur;
        end

        for i = 4:N+3
            % 选择 TV 最小的一个
            TVPol = abs(facesPol(i).right - facesPol(i+1).left) + abs(facesPol(i-1).right - facesPol(i).left);
            TVExp = abs(facesExp(i).right - facesExp(i+1).left) + abs(facesExp(i-1).right - facesExp(i).left);
            if TVPol <= TVExp
            else
                faces(i).left  = facesExp(i).left;
                faces(i).right = facesExp(i).right;
            end
        end
    end

    % 计算数值通量 (LLF)
    rhs = zeros(1, N+6);
    for i = 4:N+3
        % --------- 边界值 ----------
        % boundary(i+1/2)
        u_iphalf_L = faces(i).right;     % 来自 cell(i) 的右边界值 (^-)
        u_iphalf_R = faces(i+1).left;    % 来自 cell(i+1) 的左边界值 (^+)
    
        % boundary(i-1/2)
        u_imhalf_L = faces(i-1).right;   % 来自 cell(i-1) 的右边界值 (^-)
        u_imhalf_R = faces(i).left;      % 来自 cell(i) 的左边界值   (^+)
    
        % --------- LLF 数值通量 ----------
        % 右边界通量 (i+1/2)
        flux_iphalf = lfFlux(flux, alpha, u_iphalf_L, u_iphalf_R);

        % 左边界通量 (i-1/2)
        flux_imhalf = lfFlux(flux, alpha, u_imhalf_L, u_imhalf_R);

        % --------- 更新 RHS ----------
        rhs(i) = -(flux_iphalf - flux_imhalf) / dx;
    end



    % 截取有效区间 (对应原始 u 的 cell)
    rhs = rhs(4:N+3);
end


