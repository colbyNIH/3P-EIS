clear all; clc
syms Rsol_A Rsol_B Ra Ca Rb Cb Rs Ea Eb real positive
syms w real positive
syms s real positive

Za = Ra./(1+1i.*w.*Ra.*Ca);
Zb = subs(Za,[Ra Ca ], [Rb Cb ]);
Zs = Rs;
Zep = simplifyFraction(Rsol_A+Rsol_B+Zs.*(Za+Zb)./(Zs+Za+Zb))

Za_tf = simplifyFraction(Za./Zs.*(Rsol_A+Rsol_B+Zs-Zep))