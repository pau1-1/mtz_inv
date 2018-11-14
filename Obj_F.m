function[F, gradF]= Obj_F(sig, PSO, MTZ_new_1D, borders, Params)
d = length(sig);%dimensions
A = diag(sig);
a = 1;
eps = 1;
n = round(4+3*log(d)); %number of samples fron MC
x_rand = zeros(n, d);
[bestSolution,~,~] = PSO(MTZ_new_1D, borders, Params);
bs = MTZ_new_1D(bestSolution.xbest);
E = 0; %expected value of mean from MC samples
for i =1:n
    z = randn(1,5);
    x_rand(i, :) = bestSolution.xbest + z*A;
    E = E + MTZ_new_1D(x_rand(i, :))/n;
end;
F = log(trace(A))+ a*(max(0, E-eps*bs))^2; %objective function
dP = 2*a*max(0, E - eps*bs); %dP/dE, where P is the Penaulty function
dE = zeros(1, d); %dE/dsig, where sig is diagonal elements of matrix A
gradF = zeros(1, d); %gradient of objective function
for i= 1:d
    for j = 1:n
        x = x_rand(j,:);
        dE(i) = dE(i) - MTZ_new_1D(x_rand(j, :))*(1/2/(sig(i)+1e-8) + sig(i)*(x_rand(j, d-i+1)-bestSolution.xbest(i))^2/(prod(sig)^2+1e-8) + norm((x(end:-1:1) - bestSolution.xbest(end:-1:1))*A)^2/(sig(i)*prod(sig)^2+1e-8))/n; % dE/dsig(i)
    end;
    gradF(i) = 3/2/(sig(i)+1e-8) + dP*dE(i);
end;
gradF = -gradF;
end
