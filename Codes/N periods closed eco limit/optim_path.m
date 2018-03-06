function F = optim_path(~, BgT, Bhvec, iter, tau1, targetBg, solver)
	global chi beta R rho kappa kappa_n G alpha psi CARA id Delta CD haircut pi u h inv_h beta_borr gamma gamma_borr vars Ts Tb W Ws Wb h inv_h names Nbar Tbar horizon delev_start delev_end pi_jump tau0 allow2E privdel_end openness lambda theta
	%% Find first-best / flexible price eq'm if no default
	% Notice that Bht is B^h_{t-1} in the notes (same with Bg)
	Bh0 = Bhvec(1); BhT = Bhvec(end);

	paths = NaN(length(Bhvec),length(names));
	for j = 1:length(names)
		eval(['vars.',names{j},' = ', num2str(j),';'])
	end

	paths(:,vars.Bh) = Bhvec;
	% paths(:,vars.Bg) = Bgvec;

	% Find the no-default steady-state at time T
	% Find labor supply in the no-default steady state
	z = fsolve( @(Z) final_steady_state(Z,BhT,BgT,Tb,Ts,0,1), [Nbar,1], optimoptions('fsolve','Display','off'));
	paths(end,vars.N)		= z(1);
	paths(end,vars.Ph)		= z(2);

	% and all other variables
	z = final_steady_state(z,BhT,BgT,Tb,Ts,0,2);
	paths(end,vars.Cb)		= z(1);
	paths(end,vars.Cs)		= z(2);
	paths(end,vars.Nb)		= z(3);
	paths(end,vars.Ns)		= z(4);
	paths(end,vars.T)		= Tbar;
	paths(end,vars.q) 		= 1;
	paths(end,vars.q0)		= 1;
	paths(end,vars.pidef) 	= 0;
	paths(end,vars.Bg)		= BgT;
	% paths(end,vars.Bh)		= BhT;
	Us = u(paths(end,vars.Cs),paths(end,vars.Ns)); paths(end,vars.Vs) = Us;
	Ub = u(paths(end,vars.Cb),paths(end,vars.Nb)); paths(end,vars.Vb) = Ub;

	% Compute the steady-state at the end of the postdefault path
	z = fsolve( @(Z) final_steady_state(Z,BhT,(1-haircut)*BgT,Tb,Ts,0,1), [paths(end,vars.N), paths(end,vars.Ph)], optimoptions('fsolve','Display','off'));
	paths(end,vars.Nd)	= z(1);
	paths(end,vars.Phd)	= z(2);

	z = final_steady_state(z,BhT,(1-haircut)*BgT,Tb,Ts,0,2);
	paths(end,vars.Cbd)	= z(1);
	paths(end,vars.Csd)	= z(2);
	paths(end,vars.Nbd)	= z(3);
	paths(end,vars.Nsd)	= z(4);
	paths(end,vars.Td)	= G + (1-haircut) * (1-beta)/beta * BgT;
	paths(end,vars.qd) 	= 1;

	% Utility at the end of path after default
	Usd = u(paths(end,vars.Csd),paths(end,vars.Nsd)); paths(end,vars.Vsd) = Usd;
	Ubd = u(paths(end,vars.Cbd),paths(end,vars.Nbd)); paths(end,vars.Vbd) = Ubd;

	% Compute the steady-state at the end of the default path
	z = fsolve( @(Z) final_steady_state(Z,BhT,(1-haircut)*BgT,Tb,Ts,0,1), [paths(end,vars.N), paths(end,vars.Ph)], optimoptions('fsolve','Display','off'));
	paths(end,vars.Ndi)		= z(1);
	paths(end,vars.Phdi)	= z(2);

	z = final_steady_state(z,BhT,(1-haircut)*BgT,Tb,Ts,1,2);
	paths(end,vars.Cbdi)	= z(1);
	paths(end,vars.Csdi)	= z(2);
	paths(end,vars.Nbdi)	= z(3);
	paths(end,vars.Nsdi)	= z(4);
	paths(end,vars.Tdi)		= G + (1-haircut) * (1-beta)/beta * BgT;

	% Utility at the end of path after default
	Usdi = u(paths(end,vars.Csdi),paths(end,vars.Nsdi)); paths(end,vars.Vsdi) = Usdi;
	Ubdi = u(paths(end,vars.Cbdi),paths(end,vars.Nbdi)); paths(end,vars.Vbdi) = Ubdi;


	% "effective" wage for savers in default (calculated from wages plus profits)
	wd = (1 - chi - Delta) / (1 - chi);

	%% Now we start going back
	jtime = 1; doubleEuler = 0;
	while jtime < horizon
		jt = horizon-jtime;
		paths(jt, vars.Ph) = 1;
		paths(jt, vars.Phd) = 1;
		% Expectations

		% Remember that bh is B^h_{t-1} and bhp is B^h_t
		bhp = Bhvec(jt+1); bh = Bhvec(jt);
		bgp = paths(jt+1,vars.Bg);

		fut(:) = paths(jt+1,:);

		% Run solver for the period
		if jt > delev_end
			[curr,z] = outcomes_period(BgT,bgp,bh,bhp,fut);
			RVs = z(1); RVb = z(2); RVsi = z(3); RVbi = z(4); pid = z(5);
			eff_tau = curr(vars.T)/curr(vars.N);
		elseif jt >= delev_start && iter == 3	% deleveraging
			k = delev_start + .125 * 40;		% 40 = delev_end - delev_start for the immediate start case.
			dif = .15;
			if delev_start >= 10
				k = delev_start + (.15 + .025 * (delev_start - 10)) * 40;
				dif = .05;				% To prevent the borrowers Euler equation from not holding when you wait
			end
			if jt <= k
				t1 = tau0 + dif * (tau1 - tau0) + ((jt - delev_start)/(k - delev_start)) * (1-dif) * (tau1 - tau0);
			else
				t1 = tau1;
			end
			[curr,z] = target_taxrate(t1,bgp,bh,bhp,fut);
			RVs = z(1); RVb = z(2); RVsi = z(3); RVbi = z(4); pid = z(5);
		elseif jt >= privdel_end && iter == 2	% no deleveraging
			k = privdel_end + 12;
			if jt <= k
				dif = .75;
				t1 = tau0 + dif * (eff_tau - tau0) + (jt - privdel_end)/(k - privdel_end) * (1-dif) * (eff_tau - tau0);
				[curr,z] = target_taxrate(t1, bgp, bh, bhp, fut);
			else
				[curr,z] = outcomes_period(BgT,bgp,bh,bhp,fut);
			end
			RVs = z(1); RVb = z(2); RVsi = z(3); RVbi = z(4); pid = z(5);
		else
			[curr,z] = target_taxrate(tau0,bgp,bh,bhp,fut);
			RVs = z(1); RVb = z(2); RVsi = z(3); RVbi = z(4); pid = z(5);
		end

		for jc = 1:length(curr)
			paths(jt,jc) = curr(jc);
		end

		bg = paths(jt,vars.Bg);

		% Check that borrowers aren't on Euler
		if jt >= pi_jump
			Ubn = u(paths(jt+1,vars.Cb), paths(jt+1,vars.Nb));
			Ubd	= u(paths(jt+1,vars.Cbd), paths(jt+1,vars.Nbd));
			Ub	= u(paths(jt,vars.Cb), paths(jt,vars.Nb));
			if CD
				sdf_db = beta_borr * (paths(jt+1,vars.Cbd)/paths(jt,vars.Cb))^-1 * (Ubd/Ub)^((psi-1)/psi) * (paths(jt+1,vars.Vbd)/RVb)^(1/psi - gamma_borr);
				sdf_nb = beta_borr * (paths(jt+1,vars.Cb)/paths(jt,vars.Cb))^-1  * (Ubn/Ub)^((psi-1)/psi) * (paths(jt+1,vars.Vb)/RVb)^(1/psi - gamma_borr);
			end
			if 1/R < pid * sdf_db + (1-pid)*sdf_nb
				paths(jt,vars.E2) = 1e10;
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
				elseif ~solver
					disp(['Warning: Borrowers have interior solution at t = ',num2str(jt)])
					if iter == 3
						disp(['Govt deleveraging starts at t = ', num2str(delev_start)])
						for dsadsa = 1:1e8
							dsa = 1;
						end
						clear dsadsa
					end
				end
			else
				paths(jt,vars.E2) = -1e10;
			end
		end

		% % Now iterate on the default path ( here everything is regular expected utility --no risk )
		% paths(jt,vars.Csd) = paths(jt+1,vars.Csd);
		% paths(jt,vars.qd) = 1;
		% paths(jt,vars.Td) = G + (1-haircut)*kappa*bg - paths(jt,vars.qd) * (1-haircut) * (bgp-(1-rho)*bg);
		%
		% % Output from savers' budget constraint
		% paths(jt,vars.Nd) = paths(jt,vars.Csd) + paths(jt,vars.Td) - Ts + (chi/(1-chi))*(bhp/R-bh) + (1/(1-chi)) * paths(jt,vars.qd) * (1-haircut) * lambda * (bgp-(1-rho)*bg) - (1-haircut)*lambda*(kappa/(1-chi)) * bg;
		% paths(jt,vars.Cbd) = paths(jt,vars.Nd) + bhp/R - bh - paths(jt,vars.Td) + Tb;
		%
		% Ct = chi * paths(jt,vars.Cbd) + (1-chi) * paths(jt,vars.Csd);
		%
		% CAt = kappa * bg - paths(jt,vars.qd) * (bgp - (1-rho) * bg);
		% CAt = CAt * (1-lambda) * (1-haircut);
		%
		% if abs( paths(jt,vars.Nd) - G - Ct - CAt  ) > 1e-10
		% 	disp('Market clearing after default')
		% end
		%
		% paths(jt,vars.Nsd) = W * paths(jt,vars.Nd) / Ws;
		% paths(jt,vars.Nbd) = W * paths(jt,vars.Nd) / Wb;

		Us		= u(paths(jt,vars.Cs), paths(jt,vars.Ns));
		Usd		= u(paths(jt,vars.Csd), paths(jt,vars.Nsd));
		Usdi	= u(paths(jt,vars.Csdi), paths(jt,vars.Nsdi));
		Ub		= u(paths(jt,vars.Cb), paths(jt,vars.Nb));
		Ubd		= u(paths(jt,vars.Cbd), paths(jt,vars.Nbd));
		Ubdi	= u(paths(jt,vars.Cbdi), paths(jt,vars.Nbdi));

		if CD
			if psi == 1
				paths(jt,vars.Vs)	= exp( (1-beta) * log(Us)   + beta * log(RVs) );
				paths(jt,vars.Vsdi)	= exp( (1-beta) * log(Usdi) + beta * log(RVsi) );
				paths(jt,vars.Vsd)	= exp( (1-beta) * log(Usd)  + beta * log(paths(jt+1,vars.Vsd)) );

				paths(jt,vars.Vb)	= exp( (1-beta_borr) * log(Ub)   + beta_borr * log(RVb) );
				paths(jt,vars.Vbdi)	= exp( (1-beta_borr) * log(Ubdi) + beta_borr * log(RVbi) );
				paths(jt,vars.Vbd)	= exp( (1-beta_borr) * log(Ubd)  + beta_borr * log(paths(jt+1,vars.Vbd)) );
			else
				paths(jt,vars.Vs)	= ( (1-beta) * Us^(1-1/psi)   + beta * (RVs)^(1-1/psi) )^(1/(1-1/psi));
				paths(jt,vars.Vsdi)	= ( (1-beta) * Usdi^(1-1/psi) + beta * (RVsi)^(1-1/psi) )^(1/(1-1/psi));
				paths(jt,vars.Vsd)	= ( (1-beta) * Usd^(1-1/psi)  + beta * (paths(jt+1,vars.Vsd))^(1-1/psi) )^(1/(1-1/psi));

				paths(jt,vars.Vb)	= ( (1-beta_borr) * Ub^(1-1/psi)   + beta_borr * (RVb)^(1-1/psi) )^(1/(1-1/psi));
				paths(jt,vars.Vbdi)	= ( (1-beta_borr) * Ubdi^(1-1/psi) + beta_borr * (RVbi)^(1-1/psi) )^(1/(1-1/psi));
				paths(jt,vars.Vbd)	= ( (1-beta_borr) * Ubd^(1-1/psi)  + beta_borr * (paths(jt+1,vars.Vbd))^(1-1/psi) )^(1/(1-1/psi));
			end
		% elseif CARA
		% 	Vs(jt)	= ( (1-beta) * Us^(1-psi)  + beta * RVs^(1-psi) )^(1/(1-psi));
		% 	Vsd(jt)	= ( (1-beta) * Usd^(1-psi) + beta * Vsd(jt+1)^(1-psi) )^(1/(1-psi));
		%
		% 	Vb(jt)	= ( (1-beta) * Ub^(1-psi)  + beta_borr * RVb^(1-psi) )^(1/(1-psi));
		% 	Vbd(jt)	= ( (1-beta) * Ubd^(1-psi) + beta_borr * Vbd(jt+1)^(1-psi) )^(1/(1-psi));
		% elseif id
		% 	if psi == 1
		% 		Vs(jt) 	= exp( (1-beta) * log(Cs(jt))  + beta * log(RVs) );
		% 		Vsd(jt) = exp( (1-beta) * log(Csd(jt)) + beta * log(Vsd(jt+1)) );
		%
		% 		Vb(jt) 	= exp( (1-beta) * log(Cb(jt))  + beta_borr * log(RVb) );
		% 		Vbd(jt) = exp( (1-beta) * log(Cbd(jt)) + beta_borr * log(Vbd(jt+1)) );
		% 	else
		% 		Vs(jt) 	= ( (1-beta) * Cs(jt)^(1-psi)  + beta * RVs^(1-psi) )^(1/(1-psi));
		% 		Vsd(jt) = ( (1-beta) * Csd(jt)^(1-psi) + beta * Vsd(jt+1)^(1-psi) )^(1/(1-psi));
		%
		% 		Vb(jt) 	= ( (1-beta) * Cb(jt)^(1-psi)  + beta_borr * RVb^(1-psi) )^(1/(1-psi));
		% 		Vbd(jt) = ( (1-beta) * Cbd(jt)^(1-psi) + beta_borr * Vbd(jt+1)^(1-psi) )^(1/(1-psi));
		% 	end
		end

		% if jt > delev_end
		% 	Vs_beginning = 0;
		% 	Vsd_beginning = 0;
		% else
		% 	Vs_beginning 	= u(Cs(jt),	Ns(jt))		+ beta		* ((1-pid)*Vs_beginning + pid*Vsd_beginning);
		% 	Vsd_beginning	= u(Csd(jt),Nsd(jt))	+ beta		* Vsd_beginning;
		% end
		if ~id && paths(jt,vars.Vs) < 0
			stop
		end
		if doubleEuler
			jtime = jtime - 1;
		else
			jtime = jtime + 1;
		end
	end

	paths(:,vars.Y)   = paths(:,vars.N);
	paths(:,vars.Ydi) = paths(:,vars.Ndi) * (1-Delta);
	paths(:,vars.Yd)  = paths(:,vars.Nd);

if solver == 1
	F = targetBg - paths(pi_jump,vars.Bg);
else
	F = paths;
end
