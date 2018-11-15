function [rr] = MTZ(x, sqT)
n = length(x);
h = x(1:(round(n/2)-1));
rho = x(round(n/2):n);
%%
h = h(end:-1:1);
rho = rho(end:-1:1);
mu = 4*pi*1E-7;

lT = length(sqT);
% sqT(1:lT) = 0;
% sqT(1) = 0.01;
% for i = 2:lT
%     sqT(i) = sqT(i-1)*1.1;
% end
T = sqT.^2;
phi(1:lT) = 0;
rr(1:lT) = 0;
for i = 1:lT
    w = 2*pi/T(i);
    k1 = sqrt(-1i*w*mu/rho(1));
    
    Z(1:length(h) + 1) = 0;
    Z(1) = -1i*w*mu/k1;
    ZZ = Z(1);
    R = 1;
    for j = 1:length(h)
        k = sqrt(-1i*w*mu/rho(j+1)); %
        Z(j) = (1i*w*mu/k)*coth(-k*h(j)+acoth(ZZ*k/(1i*w*mu)));
        ZZ = Z(j);
        
        A = sqrt(rho(j)/rho(j+1)); %
        B = exp(-2*k*h(j))*(R-A)/(R+A);
        R = (1+B)/(1-B);
    end
    
    phi(i) = 180*(-atan(imag(R)/real(R)))/pi-45;
    rr(i) = (1/(w*mu))*abs(ZZ)^2; %% !!
    %i = i+1;
end;
end