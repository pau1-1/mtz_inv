clc
tic
%%
Params.NOD = 5; % number of dimensions
Params.callsMax = 1e4; % number of calls
Params.w = 0.7298; % inertia coefficient, 0.7298
Params.c1 = 1.4962; % cognitive direction multiplier, 1.4962
Params.c2 = 1.4962; % social direction multiplier, 1.4962
Params.BC_v = -0.5; % boundary conditions for v: 0 - adhesion, -1 - reflecting, ...
Params.NOP = 20; % number of particles
Params.NON = 1; % number of neighborhoods
Params.Vstart = 0; % max value of initial velocity / (bordmax-bordmin)
%%
borders.min = 10*ones(1, Params.NOD);
borders.max = 1000*ones(1, Params.NOD);
borders.Vmax = zeros(1, Params.NOD);
%%
[bestSolution, ~, ~] = PSO(@MTZ_new_1D, borders, Params);
bestFunction = MTZ_new_1D(bestSolution.xbest);
% sig0 = [1000 200 400 100 500];
sig0 = [10 20 40 10 50];
sig0 = diag(sig0);
call = 0;
step0 = 20; %gradient step
% step0 = 0.5*min(1.0/d,0.25);
sig_prev = sig0;
toc
tic
[obj_prev, grad] = Obj_F(sig_prev, bestSolution, @MTZ_new_1D, bestFunction);
toc
step = step0;
sig_new = ones(size(sig_prev));
% while step > 1e-8
% tic
% while obj > 1e-8
while norm(sig_new - sig_prev) > 1e-8
%     disp('hey')
    sig_new = sig_prev + step*grad;
%     obj_prev = obj;
    tic
    [obj, grad] = Obj_F(sig_new, bestSolution, @MTZ_new_1D, bestFunction);
    toc
    if obj > obj_prev
        sig_prev = sig_new;
        obj_prev = obj;
        call = call + 1;
        step = step0;
    else
        step = step/2;
%         obj = obj_prev;
    end;
end;
% toc
%     