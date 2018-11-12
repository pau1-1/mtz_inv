function [bestSolution, Calls, Graph] = PSO(BlackBox, rr_true)
%% Particle Swarm Optimization Algorithm
%%
Params.NOD = 3; % number of dimensions
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
%% Optimization parameters
NOD = Params.NOD; % number of dimensions
callsMax = Params.callsMax; % number of calls
w = Params.w; % inertia coefficient, 0.7298
c1 = Params.c1; % cognitive direction multiplier, 1.4962
c2 = Params.c2; % social direction multiplier, 1.4962
% c3 = Params.c3; % gradient coefficient
BC_v = Params.BC_v; % boundary conditions for v: 0 - adhesion, -1 - reflecting, ...
NOP = Params.NOP; % number of particles
NON = Params.NON; % number of neighborhoods
Vstart = Params.Vstart; % max value of initial velocity / (bordmax-bordmin)
%%
bordmin = borders.min; % max coordinate, [1 NOD]
bordmax = borders.max; % min coordinate, [1 NOD]
Vmax = borders.Vmax; % max abs value of velocity, [1 NOD]. if Vmax == 0, no velocity clamping is applied

%%
generateNeighborhood = @lbest_topology; % @lbest_topology, @ring_topology
PiN = floor(NOP/NON); % particles in neighborhood
generateX = @(d, bmin, bmax) rand([1,d]).*(bmax-bmin) + bmin;
findBestParticle = @(Swarm) min([Swarm.ObjFbest]);
% isBetter = @(a,b) a<=b;

%% Initialization
calls = 0;
t = 1; % iterations
emptyParticle.x = []; % current particle position
emptyParticle.xbest = []; % best particle position
emptyParticle.v = []; % particle velocity
emptyParticle.g = []; % antigradient direction
emptyParticle.ObjF = []; % current particle objective
emptyParticle.ObjFbest = []; % best particle objective

Swarm = repmat(emptyParticle, NOP, 1);

%-------- Initial swarm generation --------%
for i=1:NOP
    x0 = generateX(NOD, bordmin, bordmax);
    
    Swarm(i).x = x0;
    Swarm(i).xbest = x0;
    Swarm(i).v = ( Vmax + Vstart*(bordmin-bordmax).*(Vmax==0) ).*(rand(1,NOD) - 0.5)*2;
    
    [ObjF, g] = BlackBox(x0, rr_true);
    calls = calls + 1;
    
    Swarm(i).g = g;
    Swarm(i).ObjF = ObjF;
    Swarm(i).ObjFbest = Swarm(i).ObjF;
end

Calls = zeros(1, floor(callsMax/NOP));
Graph = zeros(1, floor(callsMax/NOP));

Calls(1) = calls;
[bestObjF, ~] = findBestParticle(Swarm);
Graph(1) = bestObjF;

%% Main Loop
while calls<callsMax
    t = t+1;
    
    Swarm = Swarm(randperm(NOP));
    
    %-------- velocity update --------%
    for in = 1:NON
        [SwarmN, indices, indicesN] = generateNeighborhood(Swarm,PiN,in);
        [~, bestParticle] = findBestParticle(SwarmN);
        socialDirection = SwarmN(bestParticle).xbest;
        
        for i = 1:length(indicesN)
            ii = indicesN(i);
            add1 = c1*rand([1 NOD]).*(SwarmN(ii).xbest - SwarmN(ii).x);
            add2 = c2*rand([1 NOD]).*(socialDirection - SwarmN(ii).x);
            addition = add1 + add2;
            
            SwarmN(ii).v = w*SwarmN(ii).v + addition;
            
        end
        Swarm( indices ) = SwarmN;
    end
    
    %-------- position update --------%
    for i=1:NOP
        [Swarm, calls] = updatePosition(Swarm, i, BlackBox, calls, bordmin, bordmax, BC_v);
    end

    
    %-------- графики --------%
    Calls(t) = calls;
    [bestObjF, ~] = findBestParticle(Swarm);
    Graph(t) = bestObjF;
    
end
[~, bestParticle] = findBestParticle(Swarm);
bestSolution = Swarm(bestParticle);
end

function [SwarmN, ind, indN] = lbest_topology(Swarm, PiN, i)
ind = 1+PiN*(i-1) : PiN*i;
indN = 1:length(ind);
SwarmN = Swarm( ind );
end

% function [SwarmN, D] = formRing(Swarm)
% for i=1:size(Swarm,1)
%     M(i,:) = Swarm(i).x;
% end
% i = 1;
% j = 1;
% D = 0;
% while ~all(all(isinf(M)))
%     x = M(i,:);
%     M(i,:) = Inf;
%     [i, d] = dsearchn(M,x);
%     if ~isinf(d)
%         D = D + d;
%     end
%     SwarmN(j,1) = Swarm(i);
%     j = j+1;
% end
% end
% function [SwarmN, ind, indN] = ring_topology(Swarm, PiN, i)
% ind = i-(PiN-1)/2 : i+(PiN-1)/2;
% max = size(Swarm,1);
% ind = (ind>=1 & ind<=max).*ind + (ind<1).*(ind+max) + (ind>max).*(ind-max);
% SwarmN = Swarm(ind);
% indN = ( length(ind) + 1 )/2; % there are NOP neighborhoods if topology is ring
% % velocity updates only for the one in the middle
% end