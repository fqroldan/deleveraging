function F = period_solver(cs,wd,T,Ts,Tb,new_d_h,bg,bgp,q,u,target,kind)
global chi beta rho kappa G alpha psi CARA id CD W Ws lambda

% Recover other variables implied by the guess of cs
N = cs + T - Ts + (chi/(1-chi)) * new_d_h + (1/(1-chi)) * q * lambda * (bgp-(1-rho) * bg) - (kappa/(1-chi)) * lambda * bg;

N = N / wd;

Ns = N * W / Ws;
Us = u(cs, Ns);

if kind == 1
	% Return only the difference in adjusted marginal utility
	if CD
		F = target - cs^-1 * Us^((psi-1)/psi);
	elseif CARA
		F = target - exp( - alpha * cs ) * Us^-psi;
	elseif id
		F = target - cs^(-psi);
	end
else
	% Return all the variables
	F(1) = N;
	F(2) = Ns;
	F(3) = Us;
	F(4) = N - T + (new_d_h) + Tb;	% Cb
end
