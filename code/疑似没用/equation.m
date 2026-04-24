% 清空工作区
close all;
clc;

% 全局参数定义
global g_sound1 g_sound2 g_rho10 g_rho20 g_preb;

g_sound1 = 30;
g_sound2 = 30;
g_rho10 = 1000;
g_rho20 = 1;
g_preb = 10000;

% 初始化状态变量
z1 = 0.1;
z2 = 1 - z1;

U = zeros(4, 1);
U(1) = z1 * g_rho10; 
U(2) = z2 * g_rho20;
rho = U(1) + U(2);
u = 1;
v = 1;
U(3) = rho * u;
U(4) = rho * v;

fprintf('=== 验证左右特征矩阵是否互为逆矩阵 ===\n');
fprintf('状态变量U: [%f, %f, %f, %f]\n', U(1), U(2), U(3), U(4));
fprintf('体积分数z1: %f\n', getVolumeFrac(U));
fprintf('声速: %f\n', getSoundSpeed(U));
fprintf('压力: %f\n', getPressure(U));

R = getREigenMatrix(U, nx, ny);
L = getLEigenMatrix(U, nx, ny);

L*R


% ------------- 方程形式 --------------//

function z1 = getVolumeFrac(U)
global g_sound1 g_sound2 g_rho10 g_rho20;

z1rho1 = U(1);
z2rho2 = U(2);

c1_2 = g_sound1 * g_sound1;
c2_2 = g_sound2 * g_sound2;

q = g_rho20 * c2_2 - g_rho10 * c1_2;
qs = z2rho2 * c2_2 - z1rho1 * c1_2;

r = ((q - qs) + sqrt((q - qs) * (q - qs) + 4.0 * z1rho1 * c1_2 * z2rho2 * c2_2)) / (2.0 * z2rho2 * c2_2);
z1 = r / (1.0 + r);
end

function pre = getPressure(U)
global g_sound1 g_sound2 g_rho10 g_rho20 g_preb;

z1 = getVolumeFrac(U);
z2 = 1.0 - z1;

z1rho1 = U(1);
z2rho2 = U(2);

c1_2 = g_sound1 * g_sound1;
c2_2 = g_sound2 * g_sound2;

pre1 = g_preb + c1_2 * (z1rho1 / z1 - g_rho10);
pre2 = g_preb + c2_2 * (z2rho2 / z2 - g_rho20);

% note: 过小的分母会导致数值不稳定
if z1 < z2
    pre = pre2;
else
    pre = pre1;
end
% 体积加权平均压力
% pre = z1 * pre1 + z2 * pre2;
end

function sc = getSoundSpeed(U)
global g_sound1 g_sound2;

z1 = getVolumeFrac(U);
z2 = 1.0 - z1;

z1rho1 = U(1);
z2rho2 = U(2);

rho = z1rho1 + z2rho2;
rho1 = U(1) / z1;
rho2 = U(2) / z2;

c1_2 = g_sound1 * g_sound1;
c2_2 = g_sound2 * g_sound2;

sc = sqrt(rho * (z1 / (rho1 * c1_2) + z2 / (rho2 * c2_2)));

% 另一种等效形式
% sc = sqrt(rho * (z1 * rho2 * c2_2 + z2 * rho1 * c1_2) / (rho1 * c1_2 * rho2 * c2_2));

sc = 1.0 / sc;
end

% 特征值
function lambda_max = getMaxEigenValue(U, nx, ny)
soundspeed = getSoundSpeed(U);
rho = U(1) + U(2);
u = U(3) / rho;
v = U(4) / rho;

lambda_max = abs(u * nx + v * ny) + soundspeed;
end

function Flux = getPhyFlux(Uh, nx, ny)
z1rho1 = Uh(1);
z2rho2 = Uh(2);

rho = z1rho1 + z2rho2;

u = Uh(3) / rho;
v = Uh(4) / rho;
pre = getPressure(Uh);

Flux = zeros(4, 1);

% Wakimura, H., Takagi, S., & Xiao, F. (2022).
% Symmetry-preserving enforcement of low-dissipation method based on boundary variation diminishing principle.
% Computers & Fluids, 233. doi:10.1016/j.compfluid.2021.105227
if nx > ny
    Flux(1) = z1rho1 * u;
    Flux(2) = z2rho2 * u;
    Flux(3) = rho * (u^2) + pre;
    Flux(4) = rho * u * v;
else
    Flux(1) = z1rho1 * v;
    Flux(2) = z2rho2 * v;
    Flux(3) = rho * v * u;
    Flux(4) = rho * (v^2) + pre;
end
end

% 右特征矩阵 based on {z1rho1, z2rho2, rhou, rhov}
function eigMatrix = getREigenMatrix(U, nx, ny)
global g_sound1 g_sound2;

z1rho1 = U(1);
z2rho2 = U(2);

rho = z1rho1 + z2rho2;
u = U(3) / rho;
v = U(4) / rho;

z1 = getVolumeFrac(U);
z2 = 1.0 - z1;

c = getSoundSpeed(U);

c1_2 = g_sound1 * g_sound1;
c2_2 = g_sound2 * g_sound2;

zrho1 = z1rho1 * c1_2 / (z1 * z1);
zrho2 = z2rho2 * c2_2 / (z2 * z2);

sc1_2 = c1_2 / z1 - c1_2 / z1 * zrho1 / (zrho1 + zrho2);
sc2_2 = c2_2 / z2 - c2_2 / z2 * zrho2 / (zrho1 + zrho2);

eigMatrix = zeros(4, 4);

if abs(nx) > abs(ny) % x-direction
    eigMatrix(1, 1) = z1rho1 / rho;
    eigMatrix(2, 1) = z2rho2 / rho;
    eigMatrix(3, 1) = u - c;
    eigMatrix(4, 1) = v;

    eigMatrix(1, 2) = -sc2_2 / (sc1_2 - sc2_2);
    eigMatrix(2, 2) = sc1_2 / (sc1_2 - sc2_2);
    eigMatrix(3, 2) = u;
    eigMatrix(4, 2) = v;

    eigMatrix(1, 3) = 0;
    eigMatrix(2, 3) = 0;
    eigMatrix(3, 3) = 0;
    eigMatrix(4, 3) = 1;

    eigMatrix(1, 4) = z1rho1 / rho;
    eigMatrix(2, 4) = z2rho2 / rho;
    eigMatrix(3, 4) = u + c;
    eigMatrix(4, 4) = v;
else % y-direction
    eigMatrix(1, 1) = z1rho1 / rho;
    eigMatrix(2, 1) = z2rho2 / rho;
    eigMatrix(3, 1) = u;
    eigMatrix(4, 1) = v - c;

    eigMatrix(1, 2) = -sc2_2 / (sc1_2 - sc2_2);
    eigMatrix(2, 2) = sc1_2 / (sc1_2 - sc2_2);
    eigMatrix(3, 2) = u;
    eigMatrix(4, 2) = v;

    eigMatrix(1, 3) = 0;
    eigMatrix(2, 3) = 0;
    eigMatrix(3, 3) = 1;
    eigMatrix(4, 3) = 0;

    eigMatrix(1, 4) = z1rho1 / rho;
    eigMatrix(2, 4) = z2rho2 / rho;
    eigMatrix(3, 4) = u;
    eigMatrix(4, 4) = v + c;
end
end

% 左特征矩阵 based on {z1rho1, z2rho2, rhou, rhov}
function eigMatrix = getLEigenMatrix(U, nx, ny)
global g_sound1 g_sound2;

z1rho1 = U(1);
z2rho2 = U(2);

rho = z1rho1 + z2rho2;
u = U(3) / rho;
v = U(4) / rho;

z1 = getVolumeFrac(U);
z2 = 1.0 - z1;

c = getSoundSpeed(U);
c2 = c^2;

c1_2 = g_sound1 * g_sound1;
c2_2 = g_sound2 * g_sound2;

zrho1 = z1rho1 * c1_2 / (z1 * z1);
zrho2 = z2rho2 * c2_2 / (z2 * z2);

sc1_2 = c1_2 / z1 - c1_2 / z1 * zrho1 / (zrho1 + zrho2);
sc2_2 = c2_2 / z2 - c2_2 / z2 * zrho2 / (zrho1 + zrho2);

eigMatrix = zeros(4, 4);

if abs(nx) > abs(ny) % x-direction
    eigMatrix(1, 1) = (sc1_2 + c * u) / (2.0 * c2);
    eigMatrix(1, 2) = (sc2_2 + c * u) / (2.0 * c2);
    eigMatrix(1, 3) = -1.0 / (2.0 * c);
    eigMatrix(1, 4) = 0;

    eigMatrix(2, 1) = -z2rho2 * (sc1_2 - sc2_2) / (rho * c2);
    eigMatrix(2, 2) = z1rho1 * (sc1_2 - sc2_2) / (rho * c2);
    eigMatrix(2, 3) = 0;
    eigMatrix(2, 4) = 0;

    eigMatrix(3, 1) = -v;
    eigMatrix(3, 2) = -v;
    eigMatrix(3, 3) = 0;
    eigMatrix(3, 4) = 1.0;

    eigMatrix(4, 1) = (sc1_2 - c * u) / (2.0 * c2);
    eigMatrix(4, 2) = (sc2_2 - c * u) / (2.0 * c2);
    eigMatrix(4, 3) = 1.0 / (2.0 * c);
    eigMatrix(4, 4) = 0;
else % y-direction
    eigMatrix(1, 1) = (sc1_2 + c * v) / (2.0 * c2);
    eigMatrix(1, 2) = (sc2_2 + c * v) / (2.0 * c2);
    eigMatrix(1, 3) = 0;
    eigMatrix(1, 4) = -1.0 / (2.0 * c);

    eigMatrix(2, 1) = -z2rho2 * (sc1_2 - sc2_2) / (rho * c2);
    eigMatrix(2, 2) = z1rho1 * (sc1_2 - sc2_2) / (rho * c2);
    eigMatrix(2, 3) = 0;
    eigMatrix(2, 4) = 0;

    eigMatrix(3, 1) = -u;
    eigMatrix(3, 2) = -u;
    eigMatrix(3, 3) = 1.0;
    eigMatrix(3, 4) = 0;

    eigMatrix(4, 1) = (sc1_2 - c * v) / (2.0 * c2);
    eigMatrix(4, 2) = (sc2_2 - c * v) / (2.0 * c2);
    eigMatrix(4, 3) = 0;
    eigMatrix(4, 4) = 1.0 / (2.0 * c);
end
end