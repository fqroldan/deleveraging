function F = first_best(x,RBg)
% Computes first-best steady-state allocation
global beta chi B0g G1 Gbar Nbar B1bar phi gamma_A CRRA %#ok<NUSED>
% Name variables (Notice we solve directly for kappa_n)
T = x(1); Tb = x(2); Ts = x(3); C = x(4); kappa_n = x(5); S = x(6);
R = x(7); Bh = x(8);

% Recover values
G = Gbar; N = Nbar;

% Compute conditions
F(1) = G - T + RBg * (1 - 1/R);
F(2) = chi*Tb + (1-chi)*Ts;
F(3) = C - N + Bh*(R-1)/R + T - Tb;
F(4) = C - N - (R-1)*S + T - Ts;
F(5) = (1-chi)*R*S - chi*Bh - Bg;
F(6) = R*beta - 1;
if CRRA
    F(7) = 1 - kappa_n*N^phi*C^gamma;
else
    F(7) = log(kappa_n) + phi * log(N) + gamma_A*C;
end
F(8) = Bh - B1bar*R;
end
