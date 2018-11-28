clear variables
%%
% Using PSO to find a good solution xbest
% Then starting NES algorithm to estimate the correction for xbest and sigmas
%%
x_true = [300 50 600 100 1000 300 600];
disp(['True solution: ', num2str(x_true, '%-9.1f')]);

Params.NOD = length(x_true); % number of dimensions
Params.callsMax = 1e4; % number of calls

PSOParams = Params;
PSOParams.w = 0.7298; % inertia coefficient, 0.7298
PSOParams.c1 = 1.4962; % cognitive direction multiplier, 1.4962
PSOParams.c2 = 1.4962; % social direction multiplier, 1.4962
PSOParams.BC_v = -0.5; % boundary conditions for v: 0 - adhesion, -1 - reflecting, ...
PSOParams.NOP = 20; % number of particles
PSOParams.NON = 1; % number of neighborhoods
PSOParams.Vstart = 0; % max value of initial velocity / (bordmax-bordmin)
%%
borders.min = 10*ones(1, Params.NOD);
borders.max = 1000*ones(1, Params.NOD);
borders.Vmax = zeros(1, Params.NOD);
%%
lT = 100;
sqT(1:lT) = 0;
sqT(1) = 0.01;
for i = 2:lT
    sqT(i) = sqT(i-1)*1.1;
end
rr_true = MTZ(x_true, sqT);

%%
MTZ_blackbox = @(x) MTZ_new_1D(x, rr_true, sqT);
[bestSolutionPSO, ~, ~] = PSO(MTZ_blackbox, borders, PSOParams);
xbest = bestSolutionPSO.xbest;
ObjFbest = bestSolutionPSO.ObjFbest;
disp(['PSO  solution: ', num2str(xbest, '%-9.1f')])
disp(['Objective function = ', num2str(ObjFbest)])
%%
disp('Now using this solution as mean for SNES');

SNESParams = Params;
SNESParams.mean = xbest;
SNESParams.learn_rates = [1, 3];
SNES_blackbox = @(x) max(0,abs(MTZ_blackbox(x)-ObjFbest)-1e-1);
[bestSolutionSNES, CallsSNES, GraphSNES] = SNES(SNES_blackbox, borders, SNESParams);

disp(['SNES solution: ', num2str(bestSolutionSNES.xbest, '%-9.1f')])
disp(['Final variance: ', num2str(bestSolutionSNES.var, '%-9.1f')])