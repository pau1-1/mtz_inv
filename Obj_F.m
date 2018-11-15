function[costFunction, gradF, cost2]= Obj_F(sig, bestSolution, MTZ_new_1D, bestFunction)

%% NES Parameters
d = length(sig); % dimension of target vector
A = sqrtm(sig);
% Cov = diag(sig.^2); % Covariance
Cov = sig;
logCov = log(Cov); % log Covariance for objective function

eps = 1; % difference coef
n = round(4+3*log(d)); %number of samples fron MC
%% Calculating CostFunction
xbest = bestSolution.xbest';
a = 100000 ; % max function multiplier
fit = zeros(1,n);
Z = randn(d, n);% 
X = repmat(xbest,1,n);
U = A*Z;
X = X + U;

for i = 1 : n
    fit(i) = MTZ_new_1D(X(:,i));
end

MeanFit = mean(fit);
first = trace(logCov);
second = (max(0, MeanFit-eps*bestFunction))^2;
costFunction = first - a * second; %objective function
cost2 = second;
%% Calculating gradient
dV = inv(Cov);
dP = 2*a*max(0, MeanFit - eps*bestFunction); %dP/dE, where P is the Penaulty function
M = dV;

dE = zeros(d);
for i = 1 : n
    dE = dE + fit(i) * (0.5 * M * ( X(:,i) - xbest )*( X(:,i) - xbest)' * M - 0.5 * M);
end

dE = dE / n;

gradF = dV + dP*dE;

end
