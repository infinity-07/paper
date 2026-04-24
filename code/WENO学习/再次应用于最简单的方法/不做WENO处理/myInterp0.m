function yy  = myInterp0(xx, xa, xb, uavem, uave)
    % 插值函数构造
    % 输入:
    %   xa, xb    区间端点
    %   h         区间长度
    %   uavem     平均值 (左移一格)
    %   uave      平均值 (当前格)
    %   uavep     平均值 (右移一格)  [当前没用到]
    %   lambda    指数基参数
    %
    % 输出:
    %   f_hat     最终得到的插值函数 (符号)
    %   xx, yy    在 [xa,xb] 上数值化的插值结果
    
    syms a0 a1 x
    
    % 单元中心
    xc = (xa + xb) / 2;
    h = xb - xa;
    
    % === 指数多项式2 基函数 ===
    varphi0 = 1;
    varphi1 = x-xc;
%     varphi1 = exp(lambda*(x-xc));
    % varphi2 = cosh(lambda*(x-xc));
    % varphi3 = (x-xc)*sinh(lambda*(x-xc));
    % varphi4 = (x-xc)*cosh(lambda*(x-xc));
    
    % 插值多项式 (这里只用到两个基函数)
    f = a0*varphi0 + a1*varphi1;
    
    % 构造积分条件
    eqa_m = int(f, x, xa-h, xb-h)/h == uavem;
    eqa_0 = int(f, x, xa,   xb)/h == uave;
    % eqa_p = int(f, x, xa+h, xb+h)/h == uavep; % 目前没用
    
    % 解方程
    s = solve([eqa_0, eqa_m], [a0, a1]);
    
    % 最终插值函数
    f_hat = simplify(expand(s.a0*varphi0 + s.a1*varphi1));
    
    yy = double(subs(f_hat, x, xx));
end
