clear all; clc
syms Rsol_A Rsol_B Ra Ca Rb Cb Rs Ea Eb taua taub real positive
syms R1 C1 R2 C2 tau1 tau2 real positive
syms w real positive
syms s real positive
syms gamma real positive

% SET IMPEDANCE EQUATIONS FOR EACH MEMBRANE
% za = Ra./(1+s.*taua);
za = Ra./(1+s.*Ra.*Ca);
zb = subs(za,[Ra Ca taua], [Rb Cb taub]);
zs = Rs;
z1 = subs(za,[Ra Ca taua], [R1 C1 tau1]);
z2 = subs(za,[Ra Ca taua], [R2 C2 tau2]);

% DEFINE IMPEDANCE EQUATIONS
zabs = simplify(zs.*(za+zb)./(zs+za+zb));
% zabs = simplify(subs(zabs,Ra,Rs*(gamma-Rb)));
terabs = subs(zabs,s,0);
zabs = subs(zabs,s,1i*w);
zabs_with_Rsol = simplify(Rsol_A+Rsol_B+zabs);

z12 = simplify(z1+z2);
ter12 = subs(z12,s,0);
z12_with_Rsol = simplify(Rsol_A+Rsol_B+z12);

%% calculate transfer functions
syms iapp ia ib is real 
syms va vb vs vpipette real
vrsola = iapp*Rsol_A;

ia = va/za;
vb = zb*ia;
is = (va+vb)/Rs;
va_eqn = rhs(isolate(ia==iapp-is,va));

vpipette_eqn = collect(vrsola+va_eqn,iapp)








% %% find relationships between 12 and abs
% % SEPARATE numerators and denominators
% [nabs,dabs] = numden(zabs);
% [n12,d12] = numden(z12);
% 
% nabs = collect(nabs,s);
% n12 = collect(n12,s);
% 
% dabs = collect(dabs,s);
% d12 = collect(d12,s);
% 
% % CROSS MULTIPLY
% eq1 = nabs*d12;
% eq2 = n12*dabs;
% 
% % get the coefficients of this polynomial. Subtract eq1 from eq2 to set
% % equal to 0 later
% coef_eq = coeffs(collect((eq2-eq1),s),s);
% 
% % set system of equations (soe) to 0 and simplify
% soe = simplify([0 0 0 0] == coef_eq);
% % soe(1)
% % soe(2)
% % soe(3)
% % soe(4)
% 
% % SOLVE for more complex circuit
% sol = solve(soe,[R1 tau1 R2 tau2]);
% % sol = solve(soe,[R1 C1 R2 C2]);
% 
% R1 = sol.R1(1);
% R2 = sol.R2(1);
% tau1 = sol.tau1(1);
% tau2 = sol.tau2(1);
% C1 = simplify(tau1/R1);
% C2 = simplify(tau2/R2);
% 
% % optional replace taus with RC
% R1 = simplify(subs(R1,[taua taub],[Ra*Ca Rb*Cb]));
% R2 = simplify(subs(R2,[taua taub],[Ra*Ca Rb*Cb]));
% C1 = simplify(subs(C1,[taua taub],[Ra*Ca Rb*Cb]));
% C2 = simplify(subs(C2,[taua taub],[Ra*Ca Rb*Cb]));
% tau1 = simplify(subs(tau1,[taua taub],[Ra*Ca Rb*Cb]));
% tau2 = simplify(subs(tau2,[taua taub],[Ra*Ca Rb*Cb]));
% 
% %% PLAY WITH RELATIONSHIPS
% clc
% 
% % relationships with big equations
% tau_ratio = simplify(tau1/tau2);
% tau_sum = simplify(tau1+tau2);
% tau_diff = simplify(tau1-tau2);
% tau_parallel = simplify(1/(1/tau1+1/tau2));
% 
% ter = simplify(R1+R2);
% tec = simplify(C1*C2/(C1+C2));
% R_para = simplify(R1*R2/(R1+R2));
% C_para = simplify(C1+C2);
% 
% 
% tau_mult = simplify(tau1*tau2);
% tau_mult_ter_normalize = simplify(tau_mult/ter)
% tau_mult_tec_normalize = simplify(tau_mult/tec)
% tau_mult_ter_tec_normalize = simplify(tau_mult/ter/tec)
% 
% tauEQ1 = tau_mult_ter_tec_normalize
% tauEQ2 = simplify(ter*tec)
% 
% ter_divide_tau_sum = simplify(ter/tau_sum)
% 
% % get ratio of paracellular to transcellular
% Z = tau_mult;
% dZ = simplify(Z-Z/ter);
% gamma = simplify((dZ*ter-ter*Z+Z)/(ter*Z-Z-dZ*ter))
% 
% % derive shunt to transcellular relationship
% syms x y real positive
% syms Z1 Z2 TER gamma real positive
% 
% eq1 = Z1 == x*Rs/(y+Rs);
% eq2 = Z1/TER == x/y;
% eq3 = TER == y*Rs/(y+Rs);
% eq4 = gamma == Rs/y;
% sol2 = solve([eq1 eq2 eq3 eq4],[x y Rs gamma])
% 
% 
% 
% %% try values
% clc
% ra = 1000;
% rb = 3000;
% ca = 1e-6;
% cb = 1.6e-6;
% rs = 600;
% 
% Rt = vpa(subs(ter,[Ra Rb Rs Ca Cb],[ra rb rs ca cb]),8)
% Ct = vpa(subs(tec,[Ra Rb Rs Ca Cb],[ra rb rs ca cb]),8)
% r1 = vpa(subs(R1,[Ra Rb Rs Ca Cb],[ra rb rs ca cb]),8)
% r2 = vpa(subs(R2,[Ra Rb Rs Ca Cb],[ra rb rs ca cb]),8)
% c2 = vpa(subs(C1,[Ra Rb Rs Ca Cb],[ra rb rs ca cb]),8)
% c1 = vpa(subs(C2,[Ra Rb Rs Ca Cb],[ra rb rs ca cb]),8)
% 
% t1 = vpa(subs(tau1,[Ra Rb Rs Ca Cb],[ra rb rs ca cb]),8)
% t2 = vpa(subs(tau2,[Ra Rb Rs Ca Cb],[ra rb rs ca cb]),8)
% teq1 = vpa(subs(Rt*Ct,[Ra Rb Rs Ca Cb],[ra rb rs ca cb]),8)
% teq2 = vpa(subs(tau_mult_ter_tec_normalize,[Ra Rb Rs Ca Cb],[ra rb rs ca cb]),8)
% 
% z1 = t1*t2;
% dz = z1 -z1/Rt;
% gamma1 = rs/(ra+rb)
% gamma2 = (dz*Rt-Rt*z1+z1)/(z1-dz-z1/Rt)
