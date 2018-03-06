%% Finds values for the endogenous variables assuming both agents on their Euler equation

function [ c1b,c2b,c1s,c2s,s1,n1,b1h,b1g ] = bothEuler( t1,r1,t2,guess )
global G1 B0g
[ z ] = fsolve(@(x) sys_eq(x,t1,r1,t2), guess,optimoptions('fsolve','Display','off'));
c1b = z(1); c2b = z(2); c1s = z(3); c2s = z(4); s1 = z(5); n1 = z(6); b1h = z(7);
b1g = r1 * (G1 + t1 + B0g);

end

function F = sys_eq(x,T1,R1,T2)
global beta beta_borr chi B0g B0h G1 Nbar R0S0 w1 w2 Tb Ts CRRA gamma
% Name variables
C1b = x(1); C2b = x(2); C1s = x(3); C2s = x(4); S1 = x(5); N1 = x(6); B1h = x(7);

% Compute conditions
F(1) = C2b - beta_borr*R1*C1b;
F(2) = C1b - w1*N1 - B1h/R1 + B0h + T1 - Tb;
F(3) = C2b - w2*Nbar + B1h + T2 - Tb;
F(4) = C1s - w1*N1 - R0S0 + S1 + T1 - (1-w1)/(1-chi)*N1 - Ts;
F(5) = C2s - w2*Nbar - R1*S1 + T2 - (1-w2)/(1-chi)*Nbar - Ts;
if CRRA
    F(6) = C2s - beta*R1*C1s;
else
    F(6) = C2s - C1s - log(beta*R1)/gamma;
end
F(7) = S1 - (chi/(1-chi))*B1h/R1 - (B0g-T1+G1)/(1-chi);
end
