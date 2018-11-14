clear all;
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
xbest = bestSolution.xbest;
bestFunction = MTZ_new_1D(bestSolution.xbest);
sig0 = 0.1 * eye(5);
call = 0;
step0 = 10; %gradient step
sig_prev = sig0;
[obj_best, grad_prev] = Obj_F(sig_prev, bestSolution, @MTZ_new_1D, bestFunction);
step = step0;
sig_new = zeros(size(sig_prev));
max_call = 100;
while call < max_call
    sig_new = sig_prev + step*grad_prev;
    [obj_new, grad_new] = Obj_F(sig_new, bestSolution, @MTZ_new_1D, bestFunction);
    
    if norm(sig_new - sig_prev) < 1e-8
        disp('gain is less than 1e-8')
        break
    end
    
    if obj_new > obj_best
        obj_best = obj_new;
        grad_prev = grad_new;
        sig_prev = sig_new;
        sig_best = sig_new;
        step = step0;
    else
        step = step/2;
    call = call + 1;
    end;
end;
% toc
%     