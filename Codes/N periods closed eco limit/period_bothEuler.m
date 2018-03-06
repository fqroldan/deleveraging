function F = period_bothEuler(x,w,wb,ws,T,bh,targetB,targetS,Tb,u,solver)
	global G chi CD psi beta
	% Name guess
	Cs = x(1); Cb = x(2);

	% Implied eq'm and utility
	N = chi * Cb + (1-chi) * Cs + G;

	Nb = N * w / wb;
	Ns = N * w / ws;

	Us = u(Cs, Ns);
	Ub = u(Cb, Nb);

	if solver == 1
		% Evaluate Euler eq's
		if CD
			F(1) = targetS - Cs^-1 * Us^((psi-1)/psi);
			F(2) = targetB - Cb^-1 * Ub^((psi-1)/psi);
		end
	else
		% Return all variables
		new_d = Cb - N + T - Tb;
		bhp = (new_d + bh)/beta;

		F(1) = N;
		F(2) = Ns;
		F(3) = Nb;
		F(4) = Us;
		F(5) = Ub;
		F(6) = bhp;
	end
end
