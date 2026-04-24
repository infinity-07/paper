clc,clear,close all
c1 = 1500;
c2 = 350;
rho1 = 1000;
rho2 = 1;

alpha1 = 0:1e-3:1;
alpha2 = 1-alpha1;

rho = alpha1*rho1 + alpha2*rho2;

c = ( rho .* ( alpha1./(rho1 .* c1.^2) + alpha2./(rho2 .* c2.^2) ) ).^(-0.5);

semilogy(alpha1, c)
