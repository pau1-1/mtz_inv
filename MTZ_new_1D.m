function [objF, g] = MTZ_new_1D(x, rr_true, sqT)
g = 0;
%%
rr = MTZ(x, sqT);
objF = 0;
% objF = objF + norm(abs(log(rr) - log(rr_true)));
objF = objF + sqrt((log(rr) - log(rr_true)) * (log(rr) - log(rr_true))');
% objF  =  objF + norm(abs(fi - phi_true));
end