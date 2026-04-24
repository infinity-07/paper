clc; clear; close all;

%% 符号变量
syms x xa h real
assume(h > 0);

%% ------------------ 配置部分（只需要改这里） ------------------
% 基函数（可更换为非多项式基）
varphi = [1, x, x^2];      % 例子：二次多项式基，长度自动识别为 m+1

% 模板单元相对下标（可轻松扩展）
idx = [-1, 0, 1];           % 例如：使用 u_{j-1},u_j,u_{j+1}

% 单元平均值符号变量，自动生成
u = sym("u", [1 length(idx)]);   % u(1)=u_{j-1}, u(2)=u_j, ...
u = u.';                          % 列向量

%% ---------------------------------------------------------------

% 单元端点
xa_j = @(k) xa + k*h;     % 单元 j+k 的左端点
xb_j = @(k) xa + (k+1)*h; % 右端点

m = length(varphi) - 1;   % 基函数阶数
K = length(idx);          % 模板数量

%% ------------------ 生成矩阵 A 和向量 b -----------------------
A = sym(zeros(K, m+1));

for ii = 1:K
    shift = idx(ii); 
    for jj = 1:(m+1)
        A(ii,jj) = int(varphi(jj), x, xa_j(shift), xb_j(shift));
    end
end

b = u * h;

%% ------------------ 求系数 -----------------------
c = A\b;

%% ------------------ 构造重构函数 f(x) -----------------------
f = sum(c.' .* varphi);

%% 对当前单元 (xa → xa+h) 的积分作为输出示例
result = simplify(int(f, x, xa, xa+h));

%% 输出结果
disp("A 矩阵 = ");
disp(simplify(A));

disp("b 向量 = ");
disp(simplify(b));

disp("重构函数 f(x) = ");
disp(simplify(f));

disp(simplify(subs(f,x,xb_j(0))))