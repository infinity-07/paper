clc,clear,close all
% 主要功能：
% 根据给定的模板上的原始数据，计算目标点的重构值。
% 例如：根据 u_{j-2},u_{j-1},u_{j},u_{j+1},u_{j+2} 计算 u_{j+1/2} 的值

% 这个使用的是拉格朗日函数的方法


% 定义符号变量
syms uavemm uavem uave uavep uavepp   % 不同单元的体平均值, 对应 u_{j-2}, u_{j-1}, u_{j}, u_{j+1}, u_{j+2}
syms a0 a1 a2 a3 a4                   % 插值多项式的待定系数
syms xa h x                           % 单元区间起点、单元宽度、积分变量
xb = xa+h;                            % 单元右端点
syms l2h2

% 单元中心和半宽
xc = (xa+xb)/2;       
dxc = h/2;


x_target = xb;  % 实际位置 

% % 指数多项式
syms lambda real
assume(lambda>0)
varphi0 = 1;
varphi1 = x^1;
varphi2 = x^2;
varphi3 = exp(-lambda*x);
varphi4 = exp(lambda*x);

Varphi00 = 1;
Varphi0 = int(varphi0, x);
Varphi1 = int(varphi1, x);
Varphi2 = int(varphi2, x);
Varphi3 = int(varphi3, x);
Varphi4 = int(varphi4, x);

% L_0(x_{j-3/2}) = 1; L_0(x_{j-1/2}) = 0; L_0(x_{j+1/2}) = 0; L_0(x_{j+3/2}) = 0;
% L_1(x_{j-3/2}) = 0; L_1(x_{j-1/2}) = 1; L_1(x_{j+1/2}) = 0; L_1(x_{j+3/2}) = 0;
% L_2(x_{j-3/2}) = 0; L_2(x_{j-1/2}) = 0; L_2(x_{j+1/2}) = 1; L_2(x_{j+3/2}) = 0;
% L_3(x_{j-3/2}) = 0; L_3(x_{j-1/2}) = 0; L_3(x_{j+1/2}) = 0; L_3(x_{j+3/2}) = 1;

% L_0(x) \phi_0(x_{j-3/2}-x_{j+1/2}) + L_1(x) \phi_0(x_{j-1/2}-x_{j+1/2}) + L_2(x) \phi_0(x_{j+1/2}-x_{j+1/2}) + L_3(x) \phi_0(x_{j+3/2}-x_{j+1/2}) = \phi_0(x-x_{j+1/2})
% L_1(x) \phi_1(x_{j-3/2}-x_{j+1/2}) + L_1(x) \phi_1(x_{j-1/2}-x_{j+1/2}) + L_2(x) \phi_1(x_{j+1/2}-x_{j+1/2}) + L_3(x) \phi_1(x_{j+3/2}-x_{j+1/2}) = \phi_1(x-x_{j+1/2})
% L_2(x) \phi_2(x_{j-3/2}-x_{j+1/2}) + L_1(x) \phi_2(x_{j-1/2}-x_{j+1/2}) + L_2(x) \phi_2(x_{j+1/2}-x_{j+1/2}) + L_3(x) \phi_2(x_{j+3/2}-x_{j+1/2}) = \phi_2(x-x_{j+1/2})
% L_3(x) \phi_3(x_{j-3/2}-x_{j+1/2}) + L_1(x) \phi_3(x_{j-1/2}-x_{j+1/2}) + L_2(x) \phi_3(x_{j+1/2}-x_{j+1/2}) + L_3(x) \phi_3(x_{j+3/2}-x_{j+1/2}) = \phi_3(x-x_{j+1/2})

% A = [subs(Varphi00,x,-2*h), subs(Varphi00,x,-h), subs(Varphi00,x,0), subs(Varphi00,x,h);
%     subs(Varphi0,x,-2*h),   subs(Varphi0,x,-h),   subs(Varphi0,x,0),   subs(Varphi0,x,h);
%     subs(Varphi1,x,-2*h),   subs(Varphi1,x,-h),   subs(Varphi1,x,0),   subs(Varphi1,x,h);
%     subs(Varphi2,x,-2*h),   subs(Varphi2,x,-h),   subs(Varphi2,x,0),   subs(Varphi2,x,h);];

% A = [subs(Varphi00,x,xb-3*h - x_target), subs(Varphi00,x,xb-2*h - x_target), subs(Varphi00,x,-h), subs(Varphi00,x,0), subs(Varphi00,x,h), subs(Varphi00,x,2*h);
%      subs(Varphi0,x,xb-3*h - x_target),   subs(Varphi0,x,xb-2*h - x_target),  subs(Varphi0,x,-h),   subs(Varphi0,x,0),   subs(Varphi0,x,h),   subs(Varphi0,x,2*h);
%      subs(Varphi1,x,xb-3*h - x_target),   subs(Varphi1,x,xb-2*h - x_target),  subs(Varphi1,x,-h),   subs(Varphi1,x,0),   subs(Varphi1,x,h),   subs(Varphi1,x,2*h);
%      subs(Varphi2,x,xb-3*h - x_target),   subs(Varphi2,x,xb-2*h - x_target),  subs(Varphi2,x,-h),   subs(Varphi2,x,0),   subs(Varphi2,x,h),   subs(Varphi2,x,2*h);
%      subs(Varphi3,x,xb-3*h - x_target),   subs(Varphi3,x,xb-2*h - x_target),  subs(Varphi3,x,-h),   subs(Varphi3,x,0),   subs(Varphi3,x,h),   subs(Varphi3,x,2*h);
%      subs(Varphi4,x,xb-3*h - x_target),   subs(Varphi4,x,xb-2*h - x_target),  subs(Varphi4,x,-h),   subs(Varphi4,x,0),   subs(Varphi4,x,h),   subs(Varphi4,x,2*h);];
     A = [subs(Varphi00,x,xb-3*h - xc), subs(Varphi00,x,xb-2*h - xc), subs(Varphi00,x,xb-h- xc), subs(Varphi00,x,xb - xc), subs(Varphi00,x,xb+h - xc), subs(Varphi00,x,xb+2*h - xc)
        subs(Varphi0,x,xb-3*h - xc),   subs(Varphi0,x,xb-2*h - xc),  subs(Varphi0,x,xb-h- xc),   subs(Varphi0,x,xb - xc),   subs(Varphi0,x,xb+h - xc),   subs(Varphi0,x,xb+2*h - xc)
        subs(Varphi1,x,xb-3*h - xc),   subs(Varphi1,x,xb-2*h - xc),  subs(Varphi1,x,xb-h- xc),   subs(Varphi1,x,xb - xc),   subs(Varphi1,x,xb+h - xc),   subs(Varphi1,x,xb+2*h - xc)
        subs(Varphi2,x,xb-3*h - xc),   subs(Varphi2,x,xb-2*h - xc),  subs(Varphi2,x,xb-h- xc),   subs(Varphi2,x,xb - xc),   subs(Varphi2,x,xb+h - xc),   subs(Varphi2,x,xb+2*h - xc)
        subs(Varphi3,x,xb-3*h - xc),   subs(Varphi3,x,xb-2*h - xc),  subs(Varphi3,x,xb-h- xc),   subs(Varphi3,x,xb - xc),   subs(Varphi3,x,xb+h - xc),   subs(Varphi3,x,xb+2*h - xc)
        subs(Varphi4,x,xb-3*h - xc),   subs(Varphi4,x,xb-2*h - xc),  subs(Varphi4,x,xb-h- xc),   subs(Varphi4,x,xb - xc),   subs(Varphi4,x,xb+h - xc),   subs(Varphi4,x,xb+2*h - xc)];
        


b = [subs(Varphi00, x, x-xc)
    subs(Varphi0, x, x-xc)
    subs(Varphi1, x, x-xc)
    subs(Varphi2, x, x-xc)
    subs(Varphi3, x, x-xc)
    subs(Varphi4, x, x-xc)];

L = simplify(A\b) ;

H_m3 = 0;
H_m2 = uavemm*h;
H_m1 = (uavemm+uavem)*h;
H_0  = (uavemm+uavem+uave)*h;
H_p1 = (uavemm+uavem+uave+uavep)*h;
H_p2 = (uavemm+uavem+uave+uavep+uavepp)*h;


Hm = L(1) * H_m3 + L(2) * H_m2 + L(3) * H_m1 + L(4) * H_0 + L(5) * H_p1 + L(6) * H_p2;

f_hat = simplify(expand(subs(diff(Hm,x),x,x_target)));

taylor(f_hat,h,'order',5)