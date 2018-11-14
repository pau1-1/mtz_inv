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
[bestSolution, Calls, Graph] = PSO(@MTZ_new_1D, borders, Params);
sig0 = [1000 200 400 100 500];
call = 0;
step0 = 100; %gradient step
sig_prev = sig0;
toc
tic
[obj, grad] = Obj_F(sig_prev, @PSO, @MTZ_new_1D, borders, Params);
toc
step = step0;
while step > 1e-8
    sig_new = sig_prev + step*grad;
    obj_prev = obj;
    [obj, grad] = Obj_F(sig_new, @PSO, @MTZ_new_1D, borders, Params);
    if obj<obj_prev
        sig_prev = sig_new;
        call = call + 1;
        step = step0;
    else
        step = step/2;
        obj = obj_prev;
    end;
end;
    