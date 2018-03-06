function F = final_steady_state(Z,Bht,Bgt,Tb,Ts,default_state,solver)
	n = Z(1); ph = Z(2);
	global alpha beta chi Delta kappa_n G haircut openness lambda rho kappa epsh
	R = 1/beta; r = R-1; ph0 = 1;

	if default_state
		deltadef = Delta;
	else
		deltadef = 0;
	end

	Y = n * (1-deltadef);
	w = (1-deltadef) * ph;
	P = ( (1 - openness)^(1/epsh) * ph^(1-1/epsh) + openness^(1/epsh) )^(epsh/(epsh-1));
	Pih = ph/ph0;

	% Recover other variables given guess n
	cb = ph/P * ( w/ph * n - G - r * (beta * Bht + Bgt) + Tb );
	cs = ph/P * ( w/ph * n - G + (chi/(1-chi)) * r * (beta * Bht + ((chi+lambda-1)/chi) * Bgt) + Ts );

	C = chi * cb + (1-chi) * cs;

	% nb = κ_n - (1-α)/α * cb/wb, but wb = w*N/nb
	% so nb = κ_n / ( 1 + (1-α)/α * cb/((1-Δ)*n) )
	nb = kappa_n / ( 1 + (1-alpha)/alpha * cb/(n*w) );
	wb = w * n / nb;

	% same thing for savers
	ns = kappa_n / ( 1 + (1-alpha)/alpha * cs/(n*w) );
	ws = w * n / ns;

	if solver == 1
		F(1) = nb^chi * ns^(1-chi) - n;
		F(2) = (1-lambda) * (Bgt - (1-rho) * Bgt/Pih) + Y - G - P/ph * C - kappa * (1-lambda) * Bgt/Pih;
	else
		F(1) = cb;
		F(2) = cs;
		F(3) = nb;
		F(4) = ns;
		F(5) = wb;
		F(6) = ws;
		F(7) = P;
	end
