function [F,RV] = outcomes_period(bg,bgp,bh,bhp,fut)
global pi u beta beta_borr psi gamma gamma_borr vars haircut kappa rho Ts Tb CD h inv_h R G chi W Ws Wb lambda theta Delta

pid = pi(bgp);
curr(vars.pidef) = pid;

Pi_h = 1;%fut(vars.Ph);		% expected inflation
Pi_hd= 1;%fut(vars.Phd);

Us = u(fut(vars.Cs),fut(vars.Ns));
Usd = u(fut(vars.Csd),fut(vars.Nsd));
Usdi = u(fut(vars.Csdi), fut(vars.Nsdi));
Ub = u(fut(vars.Cb),fut(vars.Nb));
Ubd = u(fut(vars.Cbd),fut(vars.Nbd));
Ubdi = u(fut(vars.Cbdi),fut(vars.Nbdi));

% Start with the no-default path
% Find discount factors
if CD
	% Stuff for the sdf
	if gamma == 1
		EhVs = pid * log(fut(vars.Vsdi)) + (1-pid) * log(fut(vars.Vs));
		RVs = exp(EhVs);
		EhVsi = theta * log(fut(vars.Vsd)) + (1-theta) * log(fut(vars.Vsdi));
		RVsi = exp(EhVsi);
	else
		EhVs = pid * fut(vars.Vsdi)^(1-gamma) + (1-pid) * fut(vars.Vs)^(1-gamma);
		RVs = EhVs^(1/(1-gamma));
		EhVsi = theta * fut(vars.Vsd)^(1-gamma) + (1-theta) * fut(vars.Vsdi)^(1-gamma);
		RVsi = EhVsi^(1/(1-gamma));
	end

	if gamma_borr == 1
		EhVb = pid * log(fut(vars.Vbdi)) + (1-pid) * log(fut(vars.Vb));
		RVb = exp(EhVb);
		EhVbi = theta * log(fut(vars.Vbd)) + (1-theta) * log(fut(vars.Vbdi));
		RVbi = exp(EhVbi);
	else
		EhVb = pid * fut(vars.Vbdi)^(1-gamma_borr) + (1-pid) * fut(vars.Vb)^(1-gamma_borr);
		RVb = EhVb^(1/(1-gamma_borr));
		EhVbi = theta * fut(vars.Vbd)^(1-gamma_borr) + (1-theta) * fut(vars.Vbdi)^(1-gamma_borr);
		RVbi = EhVbi^(1/(1-gamma_borr));
	end

	ucpd 	= fut(vars.Csd)^-1 * Usd^((psi-1)/psi);
	ucpid 	= fut(vars.Csd)^-1 * Usd^((psi-1)/psi) * (fut(vars.Vsd) / RVsi)^(1/psi - gamma);
	ucpidi 	= fut(vars.Csdi)^-1 * Usdi^((psi-1)/psi) * (fut(vars.Vsdi) / RVsi)^(1/psi - gamma);
	ucpdi 	= fut(vars.Csdi)^-1 * Usdi^((psi-1)/psi) * (fut(vars.Vsdi) / RVs)^(1/psi - gamma);
	ucp  	= fut(vars.Cs) ^-1 * Us ^((psi-1)/psi) * (fut(vars.Vs) / RVs)^(1/psi - gamma);

	ucpdb 	= fut(vars.Cbd)^-1 * Ubd^((psi-1)/psi) * (fut(vars.Vbd) / RVb)^(1/psi - gamma_borr);
	ucpbi 	= fut(vars.Cbdi)^-1 * Ubdi^((psi-1)/psi) * (fut(vars.Vbdi) / RVbi)^(1/psi - gamma_borr);
	ucpb  	= fut(vars.Cb) ^-1 * Ub ^((psi-1)/psi) * (fut(vars.Vb) / RVb)^(1/psi - gamma_borr);

	MgVu	= beta * R * (pid/Pi_hd * ucpdi  + (1-pid)/Pi_h   * ucp) ;
	MgVui	= beta * R * (theta/Pi_hd * ucpid + (1-theta)/Pi_h * ucpidi) ;
	MgVud	= beta * R * (1/Pi_hd * ucpd) ;
	MgVub	= beta_borr * R * (pid/Pi_hd * ucpbi  + (1-pid)/Pi_h   * ucpb);
	MgVubi	= beta_borr * R * (theta/Pi_hd * ucpdb + (1-theta)/Pi_h * ucpbi);

	sdf_d = beta * fut(vars.Csd)^-1 * Usd^((psi-1)/psi) * (fut(vars.Vsd) / RVs)^(1/psi - gamma) / MgVu;
	sdf_di = beta * fut(vars.Csdi)^-1 * Usdi^((psi-1)/psi) * (fut(vars.Vsdi) / RVsi)^(1/psi - gamma) / MgVu;
	sdf_n = beta * fut(vars.Cs) ^-1 * Us ^((psi-1)/psi) * (fut(vars.Vs) / RVs)^(1/psi - gamma) / MgVu;
end
% This returns MgVu, MgVub, sdf_d, sdf_n.

if 0 && ~run_outside && jt < show_horizon
	disp(['sdf_d, sdf_n = ',num2str(sdf_d),', ', num2str(sdf_n)])
	disp(['Vd/Vn = ', num2str(fut(vars.Vsd)/fut(vars.Vs))])
	disp(['Cd, Cn = ', num2str(fut(vars.Csd)),', ', num2str(fut(vars.Cs))])
end

curr(vars.q)	= pid * sdf_di * (1-haircut) * (kappa + (1-rho)) + (1-pid) * sdf_n * (kappa + (1-rho) * fut(vars.q)) ;
curr(vars.q0)	= pid * beta  * (1-haircut) * (kappa + (1-rho)) + (1-pid) * beta  * (kappa + (1-rho) * fut(vars.q0)) ;

qdi = 1;

curr(vars.T)	= G + kappa * bg - curr(vars.q) * (bgp - (1-rho)*bg);

% Guess savers consumption
if CD
	guessCs	 = fut(vars.Cs);
	guessCsdi = fut(vars.Csdi);
	guessCsd = fut(vars.Csd);
elseif CARA
	guessCs	 = fut(vars.Cs);
	guessCsdi = fut(vars.Csdi);
elseif id
	guessCs  = MgVu^(-1/psi);
	guessCsdi = MgVui^(-1/psi);
end

% use solver = 1 to return difference
if abs(period_solver(guessCs,1,curr(vars.T),Ts,Tb,beta*bhp-bh,bg,bgp,curr(vars.q),u,MgVu,1)) > 1e-12
	% disp('wrong guess')
	[curr(vars.Cs), fv, flag]	= fzero(@(c) period_solver(c,1,curr(vars.T),Ts,Tb,beta*bhp-bh,bg,bgp,curr(vars.q),u,MgVu,1), guessCs, optimset('TolFun',1e-12));
	% disp(curr(vars.Cs))
else
	% if ~run_outside
	% 	disp('right guess')
	% end
	curr(vars.Cs) = guessCs;
end

% use solver = 0 to return other variables
z = period_solver(curr(vars.Cs),1,curr(vars.T),Ts,Tb,beta*bhp-bh,bg,bgp,curr(vars.q),u,MgVu,0);
curr(vars.N)	= z(1);
curr(vars.Ns)	= z(2);
Us				= z(3);
curr(vars.Cb)	= z(4);

Ct = chi * curr(vars.Cb) + (1-chi) * curr(vars.Cs);

CAt = kappa * bg - curr(vars.q) * (bgp - (1-rho) * bg);
CAt = CAt * (1-lambda);

if abs( curr(vars.N) - G - Ct - CAt  ) > 1e-10
	disp('Market clearing')
end

curr(vars.Nb)	= curr(vars.N) * W / Wb;
curr(vars.Bg)	= bg;
curr(vars.Bh)	= bh;

% Now find the intermediate default state
% use solver = 1 to return difference

% "effective" wage for savers in default (calculated from wages plus profits)
wd = (1 - chi - Delta) / (1 - chi);

curr(vars.Tdi)	= G + (1-haircut) * kappa * bg - (1-haircut) * qdi * (bgp - (1-rho)*bg);

if abs(period_solver(guessCsdi,wd,curr(vars.Tdi),Ts,Tb,beta*bhp-bh,(1-haircut)*bg,(1-haircut)*bgp,qdi,u,MgVui,1)) > 1e-12
	% disp('wrong guess')
	[curr(vars.Csdi), fv, flag]	= fzero(@(c) period_solver(c,wd,curr(vars.Tdi),Ts,Tb,beta*bhp-bh,(1-haircut)*bg,(1-haircut)*bgp,qdi,u,MgVui,1), guessCs, optimset('TolFun',1e-12));
	% disp(curr(vars.Cs))
else
	% if ~run_outside
	% 	disp('right guess')
	% end
	curr(vars.Csdi) = guessCsdi;
end

% use solver = 0 to return other variables
z = period_solver(curr(vars.Csdi),wd,curr(vars.Tdi),Ts,Tb,beta*bhp-bh,(1-haircut)*bg,(1-haircut)*bgp,qdi,u,MgVui,0);
curr(vars.Ndi)	= z(1);
curr(vars.Nsdi)	= z(2);
Usdi			= z(3);
curr(vars.Cbdi)	= z(4);

curr(vars.Nbdi)	= curr(vars.Ndi) * W / Wb;

Cti = chi * curr(vars.Cbdi) + (1-chi) * curr(vars.Csdi);

CAti = kappa * bg - qdi * (bgp - (1-rho) * bg);
CAti = CAti * (1-lambda) * (1-haircut);

if abs( curr(vars.Ndi) * (1-Delta) - G - Cti - CAti  ) > 1e-10
	disp(curr(vars.Ndi)*(1-Delta))
	disp(G + Cti + CAti)
	disp('Market clearing after default')
end

% Finally, in the post-default state
qd = 1;
curr(vars.qd) = qd;
curr(vars.Td)	= G + (1-haircut) * kappa * bg - (1-haircut) * qd * (bgp - (1-rho)*bg);

if abs(period_solver(guessCsd,1,curr(vars.Td),Ts,Tb,beta*bhp-bh,(1-haircut)*bg,(1-haircut)*bgp,qd,u,MgVud,1)) > 1e-12
	% disp('wrong guess')
	[curr(vars.Csd), fv, flag]	= fzero(@(c) period_solver(c,1,curr(vars.Td),Ts,Tb,beta*bhp-bh,(1-haircut)*bg,(1-haircut)*bgp,qd,u,MgVud,1), guessCs, optimset('TolFun',1e-12));
	% disp(curr(vars.Cs))
else
	% if ~run_outside
	% 	disp('right guess')
	% end
	curr(vars.Csd) = guessCsd;
end

% use solver = 0 to return other variables
z = period_solver(curr(vars.Csd),1,curr(vars.Td),Ts,Tb,beta*bhp-bh,(1-haircut)*bg,(1-haircut)*bgp,qd,u,MgVud,0);
curr(vars.Nd)	= z(1);
curr(vars.Nsd)	= z(2);
Usd				= z(3);
curr(vars.Cbd)	= z(4);

curr(vars.Nbd)	= curr(vars.Nd) * W / Wb;

Ctd = chi * curr(vars.Cbd) + (1-chi) * curr(vars.Csd);

CAtd = kappa * bg - qd * (bgp - (1-rho) * bg);
CAtd = CAtd * (1-lambda) * (1-haircut);

if abs( curr(vars.Nd) - G - Ctd - CAtd ) > 1e-10
	disp(curr(vars.Nd))
	disp(G + Ctd + CAtd)
	disp(curr(vars.Nd) - G - Ctd - CAtd)
	disp('Market clearing after postdefault')
end

F = curr;
RV = [RVs RVb RVsi RVbi pid];
