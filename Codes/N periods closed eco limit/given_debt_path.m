%% Find first-best / flexible price eq'm if no default
% Notice that Bht is B^h_{t-1} in the notes (same with Bg)

% Find steady-state transfers. Cs0 = Cb0 and χTb + (1-χ)Ts = 0
Ts = - (chi/(1-chi)) * (1-beta)/beta * (Bh0/R + Bg0);
Tb = (1-beta)/beta * (Bh0/R + Bg0);

% Ts = 0; Tb = 0;

Cs0		= 1 + (chi/(1-chi)) * (1-beta)/beta * (Bh0/R + Bg0) + Ts;
Cb0 	= 1 - (1-beta)/beta * (Bh0/R + Bg0) + Tb;

Bht = Bhvec(end); Bgt = Bgvec(end);


dist_labor = @(x,cb,cs)... % x = [nb,ns,wb,ws]
				[x(1) - kappa_n + (1-alpha)/alpha * cb/x(3),...
				 x(2) - kappa_n + (1-alpha)/alpha * cs/x(4),...
				 x(1) * x(3) - x(2) * x(4),...
				 x(3)^chi * x(4)^(1-chi) - 1];

% Get initial distribution of labor and wages
% When Cs0 = Cb0 we should have Nb = Ns = N0 and Ws = Wb = W0
z = fsolve(@(x) dist_labor(x,Cb0,Cs0),[Nbar, Nbar, Wbar, Wbar], optimoptions('fsolve','Display','off'));
if flag ~= 1, disp('flag'), end
Wb = z(3); Ws = z(4);

% Prepare variables
Cs 		= NaN(size(Bhvec));
names = {'Cb', 'N', 'T', 'q', 'pidef', 'Nb','Ns','q0',...
		 'Cbd','Nd','Td','qd','Csd','Nbd','Nsd','ret_def','ret_0','ret_EZ_def','ret_EZ_0','E2', 'Vs', 'Vsd', 'Vb', 'Vbd'};
for j = 1:length(names)
	eval([names{j},' = Cs;']);
end

% Find the no-default steady-state at time T
% Find labor supply in the no-default steady state
N(end) = fsolve( @(n) final_steady_state(n,Bht,Bgt,Tb,Ts,0,1), 1, optimoptions('fsolve','Display','off'));

% and all other variables
z = final_steady_state(N(end),Bht,Bgt,Tb,Ts,0,2);
Cb(end)	= z(1);
Cs(end)	= z(2);
Nb(end)	= z(3);
Ns(end)	= z(4);
T(end)	= Tbar;
q(end) 	= 1;
q0(end)		= 1;
pidef(end) 	= 0;
Us = u(Cs(end),Ns(end)); Vs(end) = Us;
Ub = u(Cb(end),Nb(end)); Vb(end) = Ub;

% Compute the steady-state at the end of the default path
Nd(end) = fsolve( @(n) final_steady_state(n,Bht,Bgt,Tb,Ts,1,1), N(end), optimoptions('fsolve','Display','off'));

z = final_steady_state(Nd(end),Bht,Bgt,Tb,Ts,1,2);
Cbd(end)	= z(1);
Csd(end)	= z(2);
Nbd(end)	= z(3);
Nsd(end)	= z(4);
Td(end)		= G + (1-haircut) * (1-beta)/beta * Bgvec(end);
qd(end) 	= 1;


% Utility at the end of path after default
Usd = u(Csd(end),Nsd(end)); Vsd(end) = Usd;
Ubd = u(Cbd(end),Nbd(end)); Vbd(end) = Ubd;

if id
	Vs = Cs(end); Vsd = Csd(end);
	Vb = Cb(end); Vbd = Cbd(end);
end

% tilting function
h = @(v) v^(1-gamma_A);
inv_h = @(h) h^(1/(1-gamma_A));
if gamma_A == 1
	h = @(v) log(v);
	inv_h = @(h) exp(h);
end

% "effective" wage for savers in default
wd = (1 - chi - Delta) / (1 - chi);

%% Now we start going back
jtime = 1; doubleEuler = 0;
while jtime < horizon
	jt = horizon-jtime;
	% Expectations
	csd = Csd(jt+1);
	cs 	= Cs(jt+1);
	cbd = Cbd(jt+1);
	cb = Cb(jt+1);

	Us = u(cs,Ns(jt+1));
	Usd = u(csd,Nsd(jt+1));
	Ub = u(cb,Nb(jt+1));
	Ubd = u(cbd,Nbd(jt+1));

	if doubleEuler == 1
		bh = bhp;
	else
		% Remember that bh is B^h_{t-1} and bhp is B^h_t
		bh = Bhvec(jt);
	end
	bhp = Bhvec(jt+1);
	doubleEuler = 0;
	bg = Bgvec(jt); bgp = Bgvec(jt+1);
	pid = pi(bgp);
	pidef(jt) = pid;

	% Start with the no-default path
	if id
		EhVs = pid * h(Vsd(jt+1)) + (1-pid) * h(Vs(jt+1));
		EhVb = pid * h(Vbd(jt+1)) + (1-pid) * h(Vb(jt+1));

		RVs = inv_h(EhVs);
		RVb = inv_h(EhVb);

		ucpd	= csd^-psi * (Vsd(jt+1)/RVs)^(psi - gamma_A);
		ucp 	= cs ^-psi * (Vs(jt+1) /RVs)^(psi - gamma_A);

		MgVu	= pid * ucpd + (1-pid) * ucp;

		sdf_d 	= beta * csd^-psi / MgVu * (Vsd(jt+1)/RVs)^(psi - gamma_A);
		sdf_n 	= beta * cs ^-psi / MgVu * (Vs(jt+1) /RVs)^(psi - gamma_A);
	elseif CARA
		% Stuff for the sdf
		EhVs = pid * h(Vsd(jt+1)) + (1-pid) * h(Vs(jt+1));
		EhVb = pid * h(Vbd(jt+1)) + (1-pid) * h(Vb(jt+1));

		RVs = inv_h(EhVs);
		RVb = inv_h(EhVb);

		ucpd	= exp( - alpha * csd ) * (Usd)^-psi * (Vsd(jt+1) / RVs)^(psi - gamma_A);
		ucp 	= exp( - alpha * cs  ) * (Us )^-psi * (Vs(jt+1)  / RVs)^(psi - gamma_A);

		MgVu = pid * ucpd + (1-pid) * ucp;

		sdf_d = beta * exp( - alpha * csd ) * (Usd)^-psi * (Vsd(jt+1) / RVs)^(psi - gamma_A) / MgVu;
		sdf_n = beta * exp( - alpha * cs  ) * (Us )^-psi * (Vs(jt+1)  / RVs)^(psi - gamma_A) / MgVu;
	elseif CD
		% Stuff for the sdf
		EhVs = pid * h(Vsd(jt+1)) + (1-pid) * h(Vs(jt+1));
		EhVb = pid * h(Vbd(jt+1)) + (1-pid) * h(Vb(jt+1));

		RVs = inv_h(EhVs);
		RVb = inv_h(EhVb);

		ucpd 	= csd^-1 * Usd^((psi-1)/psi) * (Vsd(jt+1)/RVs)^(1/psi - gamma_A);
		ucp  	= cs ^-1 * Us ^((psi-1)/psi) * (Vs(jt+1) /RVs)^(1/psi - gamma_A);

		ucpdb 	= cbd^-1 * Ubd^((psi-1)/psi) * (Vbd(jt+1)/RVb)^(1/psi - gamma_A);
		ucpb  	= cb ^-1 * Ub ^((psi-1)/psi) * (Vb(jt+1) /RVb)^(1/psi - gamma_A);

		MgVu	= beta * R * (pid * ucpd  + (1-pid) * ucp) ;
		MgVub	= beta_borr * R * (pid * ucpdb + (1-pid) * ucpb);

		sdf_d = beta * csd^-1 * Usd^((psi-1)/psi) * (Vsd(jt+1)/RVs)^(1/psi - gamma_A) / MgVu;
		sdf_n = beta * cs ^-1 * Us ^((psi-1)/psi) * (Vs(jt+1) /RVs)^(1/psi - gamma_A) / MgVu;
	end
	if 1==0 && ~run_outside && jt < show_horizon
		disp(['sdf_d, sdf_n = ',num2str(sdf_d),', ', num2str(sdf_n)])
		disp(['Vd/Vn = ', num2str(Vsd/Vs)])
		disp(['Cd, Cn = ', num2str(csd),', ', num2str(cs)])
		if jt == 1
			disp(iter)
			if iter == 2
				disp(['CsT = ', num2str(Cs(end))])
				disp(['CsTd = ', num2str(Csd(end))])
			end
		end
	end

	q(jt)	= pid * sdf_d * (1-haircut) * (kappa + (1-rho)) + (1-pid) * sdf_n * (kappa + (1-rho) * q(jt+1)) ;
	q0(jt)	= pid * beta  * (1-haircut) * (kappa + (1-rho)) + (1-pid) * beta  * (kappa + (1-rho) * q0(jt+1)) ;

	ret_0(jt+1)	= (kappa + (1-rho) * q0(jt+1)) / q0(jt);
	ret_def(jt+1) = (1-haircut) * (kappa + 1-rho) / q0(jt);

	ret_EZ_0(jt+1) = (kappa + (1-rho) * q(jt+1))/q(jt);
	ret_EZ_def(jt+1) = (1-haircut) * (kappa + 1-rho) / q(jt);

	T(jt)	= G + kappa * bg - q(jt) * (bgp - (1-rho)*bg);

	% Guess savers consumption
	if id
		guessCs = MgVu^(-1/psi);
	elseif CARA
		guessCs	= Cs(jt+1);
	elseif CD
		guessCs	= Cs(jt+1);
	end

	% use kind = 1 to return difference
	if abs(period_solver(guessCs,T(jt),Ts,Tb,beta*bhp-bh,bg,bgp,q(jt),W,Ws,u,MgVu,1)) > 1e-12
		% disp('wrong guess')
		[Cs(jt), fv, flag]	= fzero(@(c) period_solver(c,T(jt),Ts,Tb,beta*bhp-bh,bg,bgp,q(jt),W,Ws,u,MgVu,1), guessCs , optimset('TolFun',1e-12));
		% disp(fv)
	else
		if ~run_outside
			% disp('right guess')
		end
		Cs(jt) = guessCs;
	end

	% use kind = 2 to return other variables
	z = period_solver(Cs(jt),T(jt),Ts,Tb,beta*bhp-bh,bg,bgp,q(jt),W,Ws,u,MgVu,2);
	N(jt)	= z(1);
	Ns(jt)	= z(2);
	Us		= z(3);
	Cb(jt)	= z(4);

	if abs(Cs(jt) * (1-chi) + Cb(jt) * chi + G - N(jt)) > 1e-10
		disp('Market clearing')
	end

	Nb(jt)	= N(jt) * W / Wb;
	Ubn = u(Cb(jt+1),Nb(jt+1));
	Ub = u(Cb(jt), Nb(jt));

	% Check that borrowers aren't on Euler
	if id || CD
		if id
			sdf_db = beta_borr * (Cbd(jt+1)/Cb(jt))^-psi * (Vbd(jt+1)/RVb)^(psi - gamma_A);
			sdf_nb = beta_borr * (Cb(jt+1)/Cb(jt))^-psi  * (Vb(jt+1)/RVb)^(psi - gamma_A);
		elseif CD
			sdf_db = beta_borr * (Cbd(jt+1)/Cb(jt))^-1 * (Ubd/Ub)^((psi-1)/psi) * (Vbd(jt+1)/RVb)^(1/psi - gamma_A);
			sdf_nb = beta_borr * (Cb(jt+1)/Cb(jt))^-1  * (Ubn/Ub)^((psi-1)/psi) * (Vb(jt+1)/RVb)^(1/psi - gamma_A);
		end
		if 1/R < pid * sdf_db + (1-pid)*sdf_nb
			E2(jt) = 1;
			if allow2E
				if ~run_outside
					disp(['Computing interior solution for borrowers at time ', num2str(jt)])
				end
				% Here we should have the bothEuler period solver and then instruct to go back one period, and save the chosen borrower debt to be the state tomorrow
				doubleEuler = 1;
				z = fsolve( @(x) period_bothEuler(x,W,Wb,Ws,T(jt),bh,MgVub,MgVu,Tb,u,1), [Cs(jt) Cb(jt)], optimoptions('fsolve','Display','off') );
				if z(2) ~= Cb(jt)
					E2(jt)	= 1;
				else
					E2(jt) = 0;
				end
				Cs(jt) = z(1); Cb(jt) = z(2);
				z = period_bothEuler(z,W,Wb,Ws,T(jt),bh,MgVub,MgVu,Tb,u,2);
				N(jt)	= z(1);
				Ns(jt)	= z(2);
				Nb(jt)	= z(3);
				Us		= z(4);
				Ub		= z(5);
				bhp		= z(6);
			else
				disp(['Warning: Borrowers have interior solution at t = ',num2str(jt)])
				disp(['Govt deleveraging starts at t = ', num2str(delev_start)])
			end
		end
	end

	% Now iterate on the default path ( here everything is regular expected utility --no risk )
	Csd(jt) = Csd(jt+1);
	qd(jt) = 1;
	Td(jt) = G + (1-haircut)*kappa*bg - qd(jt)*(1-haircut)*(bgp-(1-rho)*bg);

	% temp = Nd * (1-χ-Δ)/(1-χ) = Nd * "wd"
	temp = Csd(jt) + Td(jt) - Ts + (chi/(1-chi))*(bhp/R-bh) + (1/(1-chi))*qd(jt)*(1-haircut)*(bgp-(1-rho)*bg) - (1-haircut)*(kappa/(1-chi))*bg;
	Nd(jt) = temp / wd;
	Cbd(jt) = Nd(jt) + bhp/R - bh - Td(jt) + Tb;

	if abs(Csd(jt) * (1-chi) + Cbd(jt) * chi + G - Nd(jt)*(1-Delta)) > 1e-10
		disp('Market clearing')
	end

	Nsd(jt) = W * Nd(jt) / Ws;
	Nbd(jt) = W * Nd(jt) / Wb;

	Usd		= u(Csd(jt), Nsd(jt));
	Ubd		= u(Cbd(jt), Nbd(jt));

	if jt >= pi_jump
		if id
			if psi == 1
				Vs(jt) 	= exp( (1-beta) * log(Cs(jt))  + beta * log(RVs) );
				Vsd(jt) = exp( (1-beta) * log(Csd(jt)) + beta * log(Vsd(jt+1)) );

				Vb(jt) 	= exp( (1-beta) * log(Cb(jt))  + beta_borr * log(RVb) );
				Vbd(jt) = exp( (1-beta) * log(Cbd(jt)) + beta_borr * log(Vbd(jt+1)) );
			else
				Vs(jt) 	= ( (1-beta) * Cs(jt)^(1-psi)  + beta * RVs^(1-psi) )^(1/(1-psi));
				Vsd(jt) = ( (1-beta) * Csd(jt)^(1-psi) + beta * Vsd(jt+1)^(1-psi) )^(1/(1-psi));

				Vb(jt) 	= ( (1-beta) * Cb(jt)^(1-psi)  + beta_borr * RVb^(1-psi) )^(1/(1-psi));
				Vbd(jt) = ( (1-beta) * Cbd(jt)^(1-psi) + beta_borr * Vbd(jt+1)^(1-psi) )^(1/(1-psi));
			end
		elseif CARA
			Vs(jt)	= ( (1-beta) * Us^(1-psi)  + beta * RVs^(1-psi) )^(1/(1-psi));
			Vsd(jt)	= ( (1-beta) * Usd^(1-psi) + beta * Vsd(jt+1)^(1-psi) )^(1/(1-psi));

			Vb(jt)	= ( (1-beta) * Ub^(1-psi)  + beta_borr * RVb^(1-psi) )^(1/(1-psi));
			Vbd(jt)	= ( (1-beta) * Ubd^(1-psi) + beta_borr * Vbd(jt+1)^(1-psi) )^(1/(1-psi));
		elseif CD
			if psi == 1
				Vs(jt)	= exp( (1-beta) * log(Us)  + beta * log(RVs) );
				Vsd(jt)	= exp( (1-beta) * log(Usd) + beta * log(Vsd(jt+1)) );

				Vb(jt)	= exp( (1-beta) * log(Ub)  + beta_borr * log(RVb) );
				Vbd(jt)	= exp( (1-beta) * log(Ubd) + beta_borr * log(Vbd(jt+1)) );
			else
				Vs(jt)	= ( (1-beta) * Us^(1-1/psi)  + beta * (RVs)^(1-1/psi) )^(1/(1-1/psi));
				Vsd(jt)	= ( (1-beta) * Usd^(1-1/psi) + beta * (Vsd(jt+1))^(1-1/psi) )^(1/(1-1/psi));

				Vb(jt)	= ( (1-beta) * Ub^(1-1/psi)  + beta_borr * (RVb)^(1-1/psi) )^(1/(1-1/psi));
				Vbd(jt)	= ( (1-beta) * Ubd^(1-1/psi) + beta_borr * (Vbd(jt+1))^(1-1/psi) )^(1/(1-1/psi));
			end
		end
	elseif iter == 3 		% No risk
		if psi == 1
			Vs(jt)	= exp( (1-beta) * log(Us)  + beta * log(Vs(jt+1)) );
			Vsd(jt)	= exp( (1-beta) * log(Usd) + beta * log(Vsd(jt+1)) );

			Vb(jt)	= exp( (1-beta) * log(Ub)  + beta_borr * log(Vb(jt+1)) );
			Vbd(jt)	= exp( (1-beta) * log(Ubd) + beta_borr * log(Vbd(jt+1)) );
		else
			Vs(jt)	= ( (1-beta) * Us^(1-1/psi)  + beta * (Vs(jt+1))^(1-1/psi) )^(1/(1-1/psi));
			Vsd(jt)	= ( (1-beta) * Usd^(1-1/psi) + beta * (Vsd(jt+1))^(1-1/psi) )^(1/(1-1/psi));

			Vb(jt)	= ( (1-beta) * Ub^(1-1/psi)  + beta_borr * (Vb(jt+1))^(1-1/psi) )^(1/(1-1/psi));
			Vbd(jt)	= ( (1-beta) * Ubd^(1-1/psi) + beta_borr * (Vbd(jt+1))^(1-1/psi) )^(1/(1-1/psi));
		end
	end

	% if jt > delev_end
	% 	Vs_beginning = 0;
	% 	Vsd_beginning = 0;
	% else
	% 	Vs_beginning 	= u(Cs(jt),	Ns(jt))		+ beta		* ((1-pid)*Vs_beginning + pid*Vsd_beginning);
	% 	Vsd_beginning	= u(Csd(jt),Nsd(jt))	+ beta		* Vsd_beginning;
	% end
	if ~id && Vs(jt) < 0
		stop
	end
	if doubleEuler
		jtime = jtime - 1;
	else
		jtime = jtime + 1;
	end
end

Cs	= vertcat(Cs0,Cs);
Csd = vertcat(0,Csd);
Cb	= vertcat(Cb0,Cb);
Cbd = vertcat(0,Cbd);
N	= vertcat(Nbar,N);
Nd	= vertcat(0,Nd);
pidef = vertcat(0,pidef);
q	= vertcat(1,q);
qd	= vertcat(1,qd);
q0	= vertcat(1,q0);
T	= vertcat(T(1),T);
Td 	= vertcat(0,Td);

Bgvec = vertcat(Bgvec(1), Bgvec);
Bhvec = vertcat(Bhvec(1), Bhvec);

Y = N; Yd = Nd * (1-Delta);
