% symbolic asymptote solver
clear all; close all; clc

syms Ra Rb Ca Cb Rs RsolA RsolB real positive
syms w real positive

Za = Ra/(1+1i*w*Ra*Ca);
Zb = Rb/(1+1i*w*Rb*Cb);

Va = (Za*Rs+RsolA*(Rs+Za+Zb))/(Rs+Za+Zb);
Vb = (Zb*Rs+RsolB*(Rs+Za+Zb))/(Rs+Za+Zb);

MemRatio = abs(Va)/abs(Vb);
slopeMeMRatio = diff(MemRatio,w)
MemRatio_low = limit(MemRatio,w,0)
MemRatio_high = limit(MemRatio,w,inf)
%%
values = [10000 1000 1000 2.4e-6 0.1e-6 0.4 0.3];
figureVals_low = vpa(subs(MemRatio_low, ...
    [Ra Rb Rs Ca Cb RsolA RsolB], ...
    values),2)
figureVals_high = vpa(subs(MemRatio_high, ...
    [Ra Rb Rs Ca Cb RsolA RsolB], ...
    values),2)
figure_ter = vpa(subs(RsolA + RsolB + Rs*(Za+Zb)/(Rs+Za+Zb), ...
    [Ra Rb Rs Ca Cb RsolA RsolB w], ...
    [values 0]),4)