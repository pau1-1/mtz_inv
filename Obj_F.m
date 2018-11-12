function[F, gradF]= Obj_F(sig, rr_true, PSO, MTZ)
d = lenght(sig);%dimensions
A = diag(sig);
a = 1;
eps = 1;
n = 4+3*log(d); %number of samples fron MC
x_rand = zeros(1, n);
[best_solution,~,~] = PSO(MTZ, rr_true);
bs = MTZ(best_solution, rr_true);
E = 0; %expected value of mean from MC samples
for i =1:n
    z = randn(1,5);
    x_rand(i) = best_solution + A*z;
    E = E + MTZ(x_rand(i), rr_true)/n;
end;
F = log(trace(A))+ a*(max(0, E-eps*bs))^2; %objective function
dP = 2*a*max(0, E - eps*bs); %dP/dE, where P is the Penaulty function
dE = zeros(1, d); %dE/dsig, where sig is diagonal elements of matrix A
gradF = zeros(1, d); %gradient of objective function
for i= 1:d
    for j = 1:n
        dE(i) = dE(i) - MTZ(x_rand(j), rr_true)*(1/2/A(i,i)+A(i,i)*(x_rand(j, i)-x(i))^2-norm(A*(x_rand(j,:) - x))^2/A(i,i))/n; % dE/dsig(i)
    end;
    gradF(i) = 3/2/A(i,i) + dP*dE(i);
end;
end
