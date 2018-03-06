%% Code for optimal speed of deleveraging project
% Francisco Roldán, NYU. 2016-17
close all
clc
if exist('welfare_gammas_running')==1 && welfare_gammas_running == 1
	run_outside = 1;
	disp('>> Calling main code from comparative statics on risk aversion and the start of deleveraging')
else
	clear all
	run_outside = 0;
	disp('>> Running main code on its own')
end
global chi beta R rho kappa kappa_n G alpha psi CARA id Delta CD haircut pi u h inv_h beta_borr gamma gamma_borr vars Ts Tb W Ws Wb names Nbar Tbar horizon delev_start delev_end pi_jump tau0 allow2E privdel_end openness lambda epsh theta

disp('>> On the Optimal Speed of Sovereign Deleveraging')
disp('>> Closed Economy Limit, Recursive Preferences, Output cost in TFP')

annual = 0;
savegraph = 0;
makedebtgraph = 0;
plot_returns = 0;

% Preset experiments:
expnames = {'privatedelev_nodefault',...	% 1	-- Government debt but no sovereign default risk, private deleveraging
			'privatedelev_default',...		% 2	-- Default risk but no government deleveraging
			'bothdelev'};					% 3	-- Default risk and government and private deleveraging

show_experiment = [1 2 3];	% First no risk, then no delev, then delev.
plotname = 'withprivatedelev';

%% Set parameter values					% Benchmarks
beta 	= .97*annual + .995*(1-annual);	% .995	 	-- Annual real interest rate = 2%
beta_borr = .97;						% .975 		-- .9 annual
chi 	= .5;							% .65		-- Proportion of borrowers

if exist('change_chi') == 1 && change_chi ~= 0
	beta_borr = beta_borr_up;
	chi = chi_up;
	disp(['>> beta_borr = ', num2str(beta_borr), ', prop of borrowers = ', num2str(chi)])
end

openness	= 0;
lambda		= 1;

if run_outside
	if OpenEco
		lambda = .8;			% .8		-- share of debt held domestically
	else
		lambda = 1;
	end
end

disp(['>> Proportion of sovereign debt held locally = ', num2str(lambda)])

CD = 1;
id = 0;
CARA = 0;
gamma = 10;
psi = 1;								%%% 1/IES
epsh = 2;

variable_speed = 1;
if variable_speed
	disp('>> Version with sovereign deleveraging paths ending at the same time')
else
	disp('>> Version with parallel sovereign deleveraging paths')
end
allow2E = 0;

if run_outside
	if change_gamma
		gamma = g_A;						% 2			-- CARA
	end
else
	pi_jump		= 5;
	delev_start = 10;
end
gamma_borr = gamma;

phi 	= 1; 								% 1			-- Frisch
disp(['>> gamma = ',num2str(gamma),'; delev_start = ', num2str(delev_start)])

haircut	= .5;							% .25 / .50
Delta = 10/100;			% Remember to uncomment mult by Y if using DeltaoverY
theta = 1/4/10;
% DeltaoverY = 1.8;			% in percentage points  - 1.8 bench, 5 high
% DeltaoverY = DeltaoverY/100;

horizon	= 50*annual + 200*(1-annual);	% 200 		-- Number of periods

disp(['>> Flexible prices and no more default risk after period ', num2str(horizon)])

% show_horizon = 16*(1 + 3*(1-annual));
if ~run_outside
	show_horizon = 64; %64;
end
rho 	= .242*annual + .067*(1-annual);	% Fraction of sovereign debt that needs refinancing (Target is half the debt repaid in 3.5 years, so average maturity = 7 years)
% rho	= .196*annual + .0452*(1-annual); % Macaulay duration of 5 years
rho 	= .2*annual + .05 * (1-annual);
kappa 	= (1-beta)/beta + rho;

R 		= 1/beta;						% Risk-free rate is exogenously fixed at Rβ = 1
r		= R - 1;
GoverY	= .2;							% .2		-- Government expenses / GDP in steady-state
BgoverY	= 1.1 * (1 + 3*(1-annual));		% 110% GDP	-- Initial govt debt
if run_outside && exist('highdebt') == 1 && highdebt == 1
	BgoverY = initBg * (1 + 3*(1-annual));
end

Cbar	= 1;							% 1 		-- Normalization
Ybar	= Cbar / (1 - GoverY - (1-lambda)*r*BgoverY);
Nbar 	= Ybar;							% CRS
Gbar 	= GoverY * Ybar;
Bg0		= BgoverY * Ybar;


alpha = (Nbar * (80/36-1) + 1)^(-1);			% Target of 36 hrs worked a week


YearGDP = Ybar * (1 + 3*(1-annual));
Bh0		= .8/chi	* YearGDP;			% 50% GDP	-- Initial household debt

T0  = Gbar + (1-beta)/beta * Bg0;
tau0  = T0 / Ybar;

delev_end		= delev_start+10*(1+3*(1-annual));
if variable_speed
	delev_end		= 5+10*(1+3*(1-annual));
end

%% Here we start with the experiment specific variables
% names = {'Cs', 'Cb', 'N', 'T', 'q', 'pidef', 'Csd', 'Cbd', 'Nd', 'Td', 'qd','Bgvec','Bhvec','Y','Yd', 'tau0', 'tau1','Ph','Pinfh'};
names = {'Cs', 'Cb', 'N', 'T', 'q', 'pidef', 'Nb', 'Ns', 'q0', 'Cbd', 'Nd', 'Td', 'qd', 'Csd', 'Nbd', 'Nsd', 'Bh', 'Bg', 'E2', 'Vs', 'Vsd', 'Vb', 'Vbd', 'Y', 'Yd', 'Ph', 'Phd', 'Csdi', 'Cbdi', 'Ndi', 'Tdi', 'Nbdi', 'Nsdi', 'Vsdi', 'Vbdi', 'Ydi', 'Phdi'};

for j = 1:length(names)
	eval([names{j},'_mat = NaN(length(expnames),horizon+1);'])
end

iter = 1;
for experiment = show_experiment
	expname = expnames{experiment};

	% Delta 	= DeltaoverY * Ybar;

	% kappa_n	= exp(-alpha - phi * log(Ybar));	% so Wbar = 1 clears the labor market
	if CD
		kappa_n = (1-alpha)/alpha + (Ybar);	% so Wbar = 1 clears the labor market at N = 1+G when C = 1
	end
	Wbar 	= 1;

	%% Forcing variables
	G 		= Gbar;
	W 		= Wbar;
	BhT 	= .7/chi	* YearGDP;			% 40% Private deleveraging
	BgT 	= BgoverY * Ybar ;			% 110% of GDP
	if experiment == 3
		BgT		= .9 * YearGDP;				% 50% Sovereign deleveraging
	end
	%% Sequences for household and government debt

	Bhvec 	= NaN(horizon,1); Bgvec = Bhvec;

	privdel_start = 1;
	privdel_end = privdel_start+5*(1+3*(1-annual))-1;
	Bhvec(privdel_start:privdel_end)	= linspace(Bh0,BhT,5*(1+3*(1-annual)));
		Bhvec(privdel_end:end)	= BhT;
		Bhvec(1:privdel_start)	= Bh0;

	if run_outside && highdebt
		Bhvec = vertcat(Bhvec(9:end), Bhvec(end)*ones(8,1));
	end

	Bgvec(delev_start+1:delev_end) 	= linspace(Bg0,BgT,delev_end-delev_start);
		Bgvec(delev_end+1:end)	= BgT;
		Bgvec(1:delev_start)	= Bg0;

	Tbar  = Gbar + (1-beta)/beta * BgT;

	%% Default probability
	pi = @(x) min(1,max(0, (5/0.5) * (.01 * x/YearGDP + (.1-.01)*(x/YearGDP - .9)*(x/YearGDP > .9)) ));
	if run_outside && (exist('pi_crisis')==1) && pi_crisis == 0
		pi = @(x) pi(x)/3;
		pi_normal = 1;
	end
	if ~annual
		pi = @(x) (1+pi(x))^(1/20) - 1;
	else
		pi = @(x) (1+pi(x))^(1/5) - 1;
	end

	if experiment == 1
		pi = @(x) 0;
	end

	nmax = 2;
	if id
		up = @(c) 1;
		u = @(c,n) c;
	elseif CARA
		% CARA marginal utility
		up = @(c) alpha * exp(-alpha*c);
		% Period utility -- Add stuff to make this always positive
		u = @(c,n) 1 + (-1) * exp(-alpha * c) + kappa_n * (nmax^(1+phi)/(1+phi) - n^(1+phi)/(1+phi));
	elseif CD
		% Period utility
		u = @(c,n) c^alpha * (kappa_n-n)^(1-alpha);
	end

	% tilting function
	h = @(v) v^(1-gamma);
	inv_h = @(h) h^(1/(1-gamma));
	if gamma == 1
		h = @(v) log(v);
		inv_h = @(h) exp(h);
	end

	%% For plots
	PLOT.c = {[0 .447 .741],[.85 .325 .098],[0 .5 0],[.929 .6 .12],[.494 .184 .556],[.8895 .4625 .109], [.3728 .5392 .1504],[.25 .25 .25],[.7 .7 .7]};
	cnames = {'blue','red','dgreen','yellow','purple','orange','pgreen','gray','lgray'};
	for jname = 1:length(cnames)
		eval(['PLOT.',cnames{jname},' = PLOT.c{',num2str(jname),'};']);
	end

	%% Run the solver to get paths
	% Prepare variable names

	% Find steady-state transfers. Cs0 = Cb0 and χTb + (1-χ)Ts = 0
	Ts = - (chi/(1-chi)) * (1-beta)/beta * (Bh0/R + lambda * Bg0);
	Tb = (1-beta)/beta * (Bh0/R + lambda * Bg0);

	Cs0		= 1 + (chi/(1-chi)) * (1-beta)/beta * (Bh0/R + ((chi+lambda-1)/chi) * Bg0) + Ts;
	Cb0 	= 1 - (1-beta)/beta * (Bh0/R + Bg0) + Tb;

	dist_labor = @(x,cb,cs)... % x = [nb,ns,wb,ws]
					[x(1) - kappa_n + (1-alpha)/alpha * cb/x(3),...
					 x(2) - kappa_n + (1-alpha)/alpha * cs/x(4),...
					 x(1) * x(3) - x(2) * x(4),...
					 x(3)^chi * x(4)^(1-chi) - 1];

	% Get initial distribution of labor and wages
	% When Cs0 = Cb0 we should have Nb = Ns = N0 and Ws = Wb = W0
	z = fsolve(@(x) dist_labor(x,Cb0,Cs0),[Nbar, Nbar, W, W], optimoptions('fsolve','Display','off'));
	if flag ~= 1, disp('flag'), end
	Wb = z(3); Ws = z(4);

	tau1 = tau0;
	if iter == 2		% no delev
		% if run_outside && delev_start ~= 5
		% 	BgT = fzero( @(bgt) optim_path(Bg0, bgt, Bhvec, iter, tau1, targetBg, 1), BgTp);
		% else
		% 	BgT = fzero( @(bgt) optim_path(Bg0, bgt, Bhvec, iter, tau1, targetBg, 1), BgT * 1.5);
		% 	BgTp = BgT;
		% end
	elseif iter == 3	% delev
		tau1	= fzero( @(t1) optim_path(Bg0, BgT, Bhvec, iter, t1, targetBg, 1), tau0 * 1.2 );
	end
	paths = optim_path(Bg0, BgT, Bhvec, iter, tau1, 0, 0);

	if iter == 1		% no risk		-- Here we also come up with the target for debt in period pi_jump
		targetBg = paths(pi_jump,vars.Bg);
	end

	for j = 1:length(names)
		eval([names{j},' = NaN(horizon,1);'])
		eval([names{j},'(:) = paths(:,vars.', names{j},');'])
	end

	Cs0		= 1 + (chi/(1-chi)) * (1-beta)/beta * (Bh0/R + lambda * Bg0) + Ts;
	Cb0 	= 1 - (1-beta)/beta * (Bh0/R + lambda * Bg0) + Tb;

	Cs	= vertcat(Cs0,Cs);
	Csd = vertcat(0,Csd);
	Csdi = vertcat(0,Csdi);
	Cb	= vertcat(Cb0,Cb);
	Cbd = vertcat(0,Cbd);
	Cbdi = vertcat(0,Cbdi);
	N	= vertcat(Nbar,N);
	Nd	= vertcat(0,Nd);
	Ndi	= vertcat(0,Ndi);
	Nb	= vertcat(Nbar,Nb);
	Nbd	= vertcat(0,Nbd);
	Nbdi= vertcat(0,Nbdi);
	Ns	= vertcat(Nbar,Ns);
	Nsd	= vertcat(0,Nsd);
	Nsdi= vertcat(0,Nsdi);
	Y	= vertcat(Ybar,Y);
	Yd  =  vertcat(0,Yd);
	Ydi =  vertcat(0,Ydi);
	pidef = vertcat(0,pidef);
	q	= vertcat(1,q);
	qd	= vertcat(1,qd);
	q0	= vertcat(1,q0);
	T	= vertcat(T0,T);
	Td 	= vertcat(0,Td);
	Tdi	= vertcat(0,Tdi);
	E2	= vertcat(0,E2);
	Ph	= vertcat(0,Ph);
	Phd	= vertcat(0,Phd);
	Phdi= vertcat(0,Phdi);

	Bg = vertcat(Bg(1), Bg);
	Bh = vertcat(Bh(1), Bh);
	% namesV = {'Vs', 'Vsd', 'Vb', 'Vbd'};
	% for j = 1:length(namesV)
	% 	eval([namesV{j},' = ', namesV{j},'(6);']);
	% end

	for j = 1:length(names)
		if isempty(strmatch(names{j},'Vs')) && isempty(strmatch(names{j},'Vsd')) && isempty(strmatch(names{j},'Vsdi')) &&  isempty(strmatch(names{j}, 'Vb')) && isempty(strmatch(names{j},'Vbd')) && isempty(strmatch(names{j},'Vbdi'))
			eval([names{j},'_mat(iter,:) = ', names{j},';']);
		end
	end

	Vs_mat1(iter) = Vs(pi_jump);
	Vsd_mat1(iter) = Vsd(pi_jump);
	Vsdi_mat1(iter) = Vsdi(pi_jump);
	Vb_mat1(iter) = Vb(pi_jump);
	Vbd_mat1(iter) = Vbd(pi_jump);
	Vbdi_mat1(iter) = Vbdi(pi_jump);

	iter = iter+1;
end

Vs_mat = Vs_mat1;
Vsd_mat = Vsd_mat1;
Vsdi_mat = Vsdi_mat1;
Vb_mat = Vb_mat1;
Vbd_mat = Vbd_mat1;
Vbdi_mat = Vbdi_mat1;

%% Make default risk positive starting at some period
for i = 2:3
	for j = 1:length(names)
		if isempty(strmatch(names{j},'Vs')) && isempty(strmatch(names{j},'Vsd')) && isempty(strmatch(names{j},'Vsdi')) &&  isempty(strmatch(names{j}, 'Vb')) && isempty(strmatch(names{j},'Vbd')) && isempty(strmatch(names{j},'Vbdi'))
			eval([names{j},'_mat(i,1:pi_jump-1) = ', names{j},'_mat(1,1:pi_jump-1);']);
			eval([names{j},'_SS = ', names{j},'_mat(1,1);']);
		end
	end
	% keep the state at the moment of jumping
	Bg_mat(i,pi_jump) = Bg_mat(1,pi_jump);
end

% Add the pre-default values to welfare
for i = 1:3
	for jtime = 0:pi_jump-1
		jt = pi_jump - jtime;
		bgp = Bgvec(jt+1);
		pid = 0*pi(bgp);
		EhVs = pid * h(Vsdi_mat(i)) + (1-pid) * h(Vs_mat(i));
		RVs = inv_h(EhVs);
		EhVsi = theta * h(Vsd_mat(i)) + (1-theta) * h(Vsdi_mat(i));
		RVsi = inv_h(EhVsi);
		EhVb = pid * h(Vbdi_mat(i)) + (1-pid) * h(Vbd_mat(i));
		RVb = inv_h(EhVb);
		EhVbi = theta * h(Vbd_mat(i)) + (1-theta) * h(Vbdi_mat(i));
		RVbi = inv_h(EhVbi);

		Cs	 = Cs_mat(1,jt);
		Csd	 = Csd_mat(1,jt);
		Csdi = Csdi_mat(1,jt);
		Cb	 = Cb_mat(1,jt);
		Cbd	 = Cbd_mat(1,jt);
		Cbdi = Cbdi_mat(1,jt);
		Ns	 = Ns_mat(1,jt);
		Nsd	 = Nsd_mat(1,jt);
		Nsdi = Nsdi_mat(1,jt);
		Nb	 = Nb_mat(1,jt);
		Nbd	 = Nbd_mat(1,jt);
		Nbdi = Nbdi_mat(1,jt);

		if CD
			if psi == 1
				Vs_mat(i)	= exp( (1-beta) * log(u(Cs,Ns))     + beta * log(RVs) );
				Vsd_mat(i)	= exp( (1-beta) * log(u(Csd,Nsd))   + beta * log(Vsd_mat(i)) );
				Vsdi_mat(i)	= exp( (1-beta) * log(u(Csdi,Nsdi)) + beta * log(RVsi) );

				Vb_mat(i)	= exp( (1-beta_borr) * log(u(Cb,Nb))  + beta_borr * log(RVb) );
				Vbd_mat(i)	= exp( (1-beta_borr) * log(u(Cbd,Nbd)) + beta_borr * log(Vbd_mat(i)) );
				Vbdi_mat(i)	= exp( (1-beta_borr) * log(u(Cbdi,Nbdi)) + beta_borr * log(RVbi) );

			else
				Vs_mat(i)	= ( (1-beta) * u(Cs,Ns)^(1-1/psi)   + beta * (RVs)^(1-1/psi) )^(1/(1-1/psi));
				Vsd_mat(i)	= ( (1-beta) * u(Csd,Nsd)^(1-1/psi) + beta * (Vsd_mat(i))^(1-1/psi) )^(1/(1-1/psi));
				Vsdi_mat(i)	= ( (1-beta) * u(Csdi,Nsdi)^(1-1/psi)   + beta * (RVsi)^(1-1/psi) )^(1/(1-1/psi));

				Vb_mat(i)	= ( (1-beta_borr) * u(Cb,Nb)^(1-1/psi)  + beta_borr * (RVb)^(1-1/psi) )^(1/(1-1/psi));
				Vbd_mat(i)	= ( (1-beta_borr) * u(Cbd,Nbd)^(1-1/psi) + beta_borr * (Vbd_mat(i))^(1-1/psi) )^(1/(1-1/psi));
				Vbdi_mat(i)	= ( (1-beta_borr) * u(Cbdi,Nbdi)^(1-1/psi)   + beta_borr * (RVbi)^(1-1/psi) )^(1/(1-1/psi));
			end
		elseif id
			Vs_mat(i)	= exp( (1-beta) * log(Cs)  + beta * log(RVs) );
			Vsd_mat(i)	= exp( (1-beta) * log(Csd) + beta * log(Vsd_mat(i)) );

			Vb_mat(i)	= exp( (1-beta_borr) * log(Cb)  + beta_borr * log(RVb) );
			Vbd_mat(i)	= exp( (1-beta_borr) * log(Cbd) + beta_borr * log(Vbd_mat(i)) );
		end
	end
end

% Welfare equivalents
for j = 1:3
	C_equiv_s(j) = Vs_mat(j)^(1/alpha) * (kappa_n - Nbar)^((alpha-1)/alpha);
	C_equiv_b(j) = Vb_mat(j)^(1/alpha) * (kappa_n - Nbar)^((alpha-1)/alpha);
end

% Output losses comparing pub del to def risk
Y_lost_beta_delev = 1/(1-beta_borr) * (N_mat(3,201) - N_mat(1,201)) / Y_SS;
Y_lost_beta_risk =  1/(1-beta_borr) * (N_mat(2,201) - N_mat(1,201)) / Y_SS;
Y_lost_beta_delev = 0; Y_lost_beta_risk = 0;
jT = 64;
for jtime = 1:jT
	jt = jT - jtime + 1;
	Y_lost_beta_delev = beta_borr * Y_lost_beta_delev + (N_mat(3,jt) - N_mat(1,jt)) / Y_SS;
	Y_lost_beta_risk = beta_borr * Y_lost_beta_risk + (N_mat(2,jt) - N_mat(1,jt)) / Y_SS;
end

%% WARNING: Everything later works because the no-default simulation is the last one

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot results
if ~run_outside
	plot_results

	if savegraph
		% eval(['print -dpng ../../Graphs/Nperiods/',plotname,'_h',num2str(100*haircut),'_D',num2str(100*DeltaoverY),'.png -r300'])
		% eval(['print -depsc ../../Graphs/Nperiods/',plotname,'_h',num2str(100*haircut),'_D',num2str(100*DeltaoverY),'.eps'])
	end

% disp(['	&def risk		', '&pub delev		', '&no def'])
% disp(['$V^s_0(0)$  &',num2str(Vs_mat(1)),'  &',num2str(Vs_mat(2)),'  &',num2str(Vs_mat(3)),'  \\'])
% disp(['$\bar{C}^s$  &',num2str(C_equiv_s(1)),'  &',num2str(C_equiv_s(2)),'  &',num2str(C_equiv_s(3)),'  \\'])
% disp(['$V^s_0(0)$  &',num2str(Vb_mat(1)),'  &',num2str(Vb_mat(2)),'  &',num2str(Vb_mat(3)),'  \\'])
% disp(['$\bar{C}^b$  &',num2str(C_equiv_b(1)),'  &',num2str(C_equiv_b(2)),'  &',num2str(C_equiv_b(3)),'  \\ '])
end
