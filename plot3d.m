SigmaL = diag(sig_best);
SigmaL = abs(real(SigmaL));

xmin = min(xbest);
xmax = max(xbest);
sigmax = max(SigmaL);
x = [xmin-sigmax:1:xmax+sigmax];


figure;
for i=1:length(xbest)
    norm = normpdf(x,xbest(i),SigmaL(i));
    plot(x,norm);
    hold on;
end