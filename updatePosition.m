function [Swarm, calls] = updatePosition(Swarm, i, BlackBox, calls, bordmin, bordmax, BC_v, rr_true, sqT)
x0 = Swarm(i).x + Swarm(i).v; % x(t+1) = x(t) + v(t+1)

% boundary check
Swarm(i).v = (x0<=bordmax & x0>=bordmin).*Swarm(i).v + BC_v.*(x0>bordmax | x0<bordmin).*Swarm(i).v;
x0 = (x0<=bordmax & x0>=bordmin).*x0 + (x0>bordmax).*bordmax + (x0<bordmin).*bordmin;

[ObjF, g] = BlackBox(x0, rr_true, sqT);
calls = calls + 1;

Swarm(i).x = x0;
Swarm(i).g = g;
Swarm(i).ObjF = ObjF;

if ObjF <= Swarm(i).ObjFbest % isBetter(ObjF, Swarm(i).ObjFbest)
    Swarm(i).ObjFbest = ObjF;
    Swarm(i).xbest = x0;
end

end