function [bestSolution, Calls, Graph] = SNES(BlackBox, borders, Params)
%%
NOD = Params.NOD;
callsMax = Params.callsMax;

bordmin = borders.min;
bordmax = borders.max;

if isfield(Params, 'NOP')
    NOP = Params.NOP;
else
    NOP = 4 + floor(3 * log(NOD));
end
if isfield(Params, 'learn_rates')
    learn_rates = Params.learn_rates;
else
    learn_rates = [1, ( 3 + log(NOD) )/( 5*sqrt(NOD) )];
end
if isfield(Params, 'mean')
    mean = Params.mean;
else
    mean = (bordmax - bordmin).*( rand(1,NOD) - 0.5 );
end
if isfield(Params, 'var')
    var = Params.var;
else
    var = 0.01*(bordmax - bordmin).*ones(1,NOD);
end
%%
calls = 0;
Calls = [];
Graph = [];
fitness = zeros(1,NOP);
%%
bestSolution.x = [];
bestSolution.ObjF = +Inf;
bestSolution.var = var;
bestSolution.xbest = [];
bestSolution.ObjFbest = +Inf;
%%
% threshold = floor(NOP/2);
% step = 1 / threshold;
% U = zeros(1, NOP);
% U(end-threshold+1:end) = step:step:1;
% U = U / sum(U);

% U = linspace(0, 2/NOP, NOP);
%%
while calls < callsMax
    S = randn(NOP,NOD); % N(0,1)
    particles = repmat(mean,NOP,1) + S.*repmat(var,NOP,1); % N(mean_vec, var_vec^2)
    
    for i=1:NOP
        fitness(i) = BlackBox(particles(i,:));
        calls = calls + 1;
    end
    
    [~, idx] = sort(fitness, 'descend');
    S = S(idx,:);
    
    U = 1 ./ (0.1*ones(1, NOP) + fitness(idx));
    U = U / sum(U);
    
    mean_grad = U * S;
    var_grad = U * (S.^2 - 1);
    
    mean = mean + learn_rates(1) * var .* mean_grad;
    var = var .* exp(learn_rates(2) / 2 * var_grad);
    
    mean = (mean<=bordmax & mean>=bordmin).*mean + (mean>bordmax).*bordmax + (mean<bordmin).*bordmin;
    
    mean_fitness = BlackBox(mean);
    calls = calls + 1;
    Graph = cat(1, Graph, mean_fitness);
    Calls = cat(1, Calls, calls);
    
    if mean_fitness < bestSolution.ObjF
        bestSolution.xbest = mean;
        bestSolution.ObjFbest = mean_fitness;
    end
    bestSolution.x = mean;
    bestSolution.ObjF = mean_fitness;
    bestSolution.var = var;
    %     disp(['variance: ', num2str(var, '%-9.1f')]);
end

end