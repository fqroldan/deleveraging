%% Stores the parameters for the Speed of Deleveraging project codes.
global beta beta_borr chi gamma_A phi kappa_n B0g B0h B1bar G1 Gbar Nbar W1 W2 Tb Ts R0S0 low_wage priv_del%#ok<NUSED>

beta 	= .97^5;						% .97^5 	-- Discount
beta_borr = .9^5;						% .5^5 		-- Higher discount
chi 	= .5;							% .5		-- Proportion of borrowers
gamma_A = 1;							% 1			-- CARA
phi 	= 1; 							% 1			-- Frisch

GoverY	= .2;							% .2		-- Government expenses / GDP in steady-state
Cbar	= 1;							% 1 		-- Normalization
Ybar	= Cbar/(1-GoverY);
Nbar 	= Ybar;							% CRS
Gbar 	= GoverY*Ybar;

B0h	    = 1 * Ybar/5;					% 50% of yearly GDP		/
B0g 	= 1.2 * Ybar/5;					% 120% of yearly GDP	/

R0S0 = B0g/(1-chi) + (chi/(1-chi)) * B0h;

kappa_n	= exp(-gamma_A - phi * log(1+Gbar));	% so Wbar = 1 clears the labor market
Wbar 	= 1;

%% Forcing variables
G1 		= Gbar;
W1 		= (1 - .1*low_wage) * Wbar;
B1bar	= (1/(1+beta) - .2*priv_del) * B0h;
% B1bar	= (1 - .65*priv_del) * B0h;					% 80% Private deleveraging
% W1 		= (1 - 0*low_wage) * Wbar;
% B1bar	= (1/(1+beta) - 0*priv_del) * B0h;

%% Specification of the pi function and utility

CRRA = 0;                                                   % Set to 1 for CRRA, 0 for CARA
if CRRA, disp('CRRA case'); end %#ok<*UNRCH>
if ~CRRA, CARA = 1; disp('CARA case'); end
%% For plots
% red (1), dark green (2), blue (3), yellow (4), purple (5), ¿pear? green (6), orange (7)
PLOT.c = {'[.85 .325 .098]','[0 .5 0]','[0 .447 .741]','[.929 .6 .12]','[.494 .184 .556]','[.3728 .5392 .1504]',...
    '[.8895 .4625 .109]'};
PLOT.red = PLOT.c{1};
PLOT.dgreen = PLOT.c{2};
PLOT.blue = PLOT.c{3};
PLOT.yellow = PLOT.c{4};
PLOT.purple = PLOT.c{5};
PLOT.pgreen = PLOT.c{6};
PLOT.orange = PLOT.c{7};
