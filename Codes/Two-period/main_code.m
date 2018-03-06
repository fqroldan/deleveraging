%% Code for the optimal speed of deleveraging project
% Francisco Roldán, NYU. 10/2016
close all; clear all; clc %#ok<CLSCR>
disp('>> On the Optimal Speed of Deleveraging')

savegraph = 0;
savediag = 0;
%% Parameters
global beta beta_borr chi gamma_A phi kappa_n B0g B0h B1bar G1 G2 Nbar W1 W2 R0S0 CRRA CARA low_wage priv_del %#ok<NUSED>

low_wage = 0;
plotname = 'no_default_low_wage';

if low_wage == 0
	priv_del = 1;
	plotname = 'no_default_priv_del';
else priv_del = 0;
end

Delev_parameters

%% Options
showdiagram = 1;

%% Set κ_n so that W = 1 is the flexible price equilibrium
kappa_n = exp(-gamma_A * Cbar)/(Nbar)^phi;


%% Utility function. CARA case.
U = @(c1,c2,n1,n2,b) (-1/gamma_A)*exp(-gamma_A*c1)-(kappa_n)*n1^(1+phi)/(1+phi)+...
	b * ((-1/gamma_A)*exp(-gamma_A*c2) - (kappa_n)*n2^(1+phi)/(1+phi));
% Marginal Utility
up = @(c) exp(-gamma_A * c);

% Figure out the distribution of labor and wages in the initial steady-state
dist_labor = @(x,ntot,cb,cs) [ntot - x(1)^chi * x(2)^(1-chi),...
							 kappa_n * x(1)^phi - x(3)*up(cb),...
							 kappa_n * x(2)^phi - x(4)*up(cs),...
							 x(1)*x(3) - x(2)*x(4)];
%
Cs0		= 1 + (chi/(1-chi)) * (1-beta)/beta * (beta * B0h + B0g);
Cb0 	= 1 - (1-beta)/beta * (beta * B0h + B0g);

z = fsolve(@(x) dist_labor(x,Nbar,Cb0,Cs0),[Nbar, Nbar, Wbar, Wbar], optimoptions('fsolve','Display','off'));
if flag ~= 1, disp('flag'), end
Nb0 = z(1); Ns0 = z(2);
n2b = Nb0; n2s = Ns0;	% Guesses for later
Wb0 = z(3); Ws0 = z(4);
w2b = Wb0; w2s = Ws0;

%% Equilibrium for each B1g
R1 = 1/beta;
% No default risk and no intermediation frictions
R1g = R1; R1h = R1;
% Choice variable is B1g
B1g = linspace(0, B0g * R1, 1001);
% B1g = linspace(0,R1*(B0g - B1bar/R1 + B0h) + G1*(1-chi)/chi,5001);

%% Prepare variables
C1s = NaN(length(B1g),1); C1b = C1s; C2s = C1s; C2b = C1s;
N1 = C1s; S1 = C1s; Sr = C1s; C1r = C1s;
T2 = C1s;
Us = C1s; Ub = C1s;
T1show = C1s;


%% Loop
eulerstop = 0; mktstop = 0;
for jb = 1:length(B1g)
	b1g = B1g(jb);
	t1 = G1 + B0g - b1g/R1g;

	c2s = 1 + (chi/(1-chi)) * (B1bar + b1g);
	% βR = 1 and no risk, so savers' consumption is constant
	c1s = c2s;

	P = 1 + (chi/(1-chi)) * (B1bar + b1g);

	% Take N so that P = F(b1g,N). In general we would define F as a function of n and then do fzero n -> F(b1g,n) - P
	temp = P + G1 + (chi/(1-chi)) * (b1g/R1g - B0g + B1bar/R1h - B0h);
	n1 = temp / (W1 + (1-W1)/(1-chi));

	n1b	= n1 * Wbar / Wb0;
	n1s	= n1 * Wbar / Ws0;
	% temp = (-1/gamma_A) * (log(beta * R1) + (-gamma_A * c2s));
	% if abs(temp-c1s) > 1e-10
	% 	disp('Warning: Euler equation fails')
	% end

    c2b   = Nbar/chi - (1-chi)/chi * c2s - Gbar/chi;
	z = fsolve(@(x) dist_labor(x,Nbar,c2b,c2s),[n2b, n2s, w2b, w2s], optimoptions('fsolve','Display','off'));
	n2b = z(1); n2s = z(2);
	w2b = z(3); w2s = z(4);

	% Check for errors
	if c2b < 0
        disp('ERROR: Negative C_2^b')
    end

	b1h = B1bar;
    c1b = W1*n1 + b1h/R1h - B0h - t1;
	if c1b < 0
		disp('ERROR: Negative C_1^b')
	end
	% What savings should be like
	s1 = (1/(1-chi)) * (b1g/R1g + chi * b1h/R1h);
    dispy = ((1-chi*W1)/(1-chi))*n1 + R0S0 - t1;
    sr = s1/dispy;
    c1r = c1s/dispy;
    temp = ((1-chi*W1)/(1-chi))*n1 + R0S0 - t1 - c1s;

	% Check that market-clearing savings satisfy the savers' budget constraint
	if abs(temp-s1) > 10^-2
		disp('market clearing')
	end

    if exp(-gamma_A*c1b) < beta_borr * (exp(-gamma_A*c2b))
        disp('Euler: Borrowers have interior solution')
    end
    t2 = b1g + Gbar;

	% Store
	C1b(jb)		= c1b;
    C1s(jb)		= c1s;
    C2b(jb)		= c2b;
    C2s(jb)		= c2s;
    S1(jb)		= s1;
    Sr(jb)		= sr;
    C1r(jb)		= c1r;
    N1(jb)		= n1;
    T1show(jb)	= t1;
    T2(jb)		= t2;

	% Welfare
	Us(jb) = U(c1s,c2s,n1s,n2s,beta);
    Ub(jb) = U(c1b,c2b,n1b,n2b,beta_borr);
end

% Find the policy that implements B1g = B0g
jdebt = find(B1g >= B0g, 1);
jdebt2 = find(B1g >= B0g*0.8,1);
jneutral = find(B1g >= B0g/(1+beta),1);

% Normalize everything by steady-state GDP
N2 = Nbar * ones(size(B1g));
normalize_vars = @(x) x / Nbar;
names = {'C1b','C1s','C2b','C2s','S1','N1','N2','B1g'};
for jname = 1:length(names)
	name = names{jname};
	eval([name, '_show = normalize_vars(',name,');'])
end

% The cross diagram
Npic = N1(jdebt);
Npic2 = N1(jdebt2);
PresentValueCurve = 1 + (chi/(1-chi)) * (B1g(:) + B1bar);
FundingCurve = (W1 + (1-W1)/(1-chi)) * Npic - G1 - (chi/(1-chi)) * (B1g(:)/R1g - B0g + B1bar/R1h - B0h);
FundingCurve2 = (W1 + (1-W1)/(1-chi)) * Npic2 - G1 - (chi/(1-chi)) * (B1g(:)/R1g - B0g + B1bar/R1h - B0h);
figure('position',[700   500   800   450]), hold on
plot(B1g_show, FundingCurve,'w-.','linewidth',1)
p2 = plot(B1g_show, FundingCurve2,'-.','linewidth',1,'color',PLOT.red)
p1 = plot(B1g_show, PresentValueCurve,'linewidth',1.5, 'color', PLOT.blue)
yl = get(gca,'ylim');
plot(B1g_show(jdebt2),yl(1),'kx','linewidth',1.5)
plot(B1g_show(jdebt2)*ones(2,1),[yl(1),PresentValueCurve(jdebt2)],'k:','linewidth',1)
l = legend([p1 p2], '$\mathcal{P}(B_1^g)$','$\mathcal{F}(B_1^g,N_1)$');
set(l,'Interpreter','Latex','location','south','orientation','horizontal');
xlabel('$B_1^g/\bar{Y}$','Interpreter','Latex')
ylabel('$C_1^s$','Interpreter','Latex')
if savediag
	eval(['print -dpng ../../Graphs/2periods/diagram_nodefault1.png -r300'])
	eval(['print -depsc ../../Graphs/2periods/diagram_nodefault1.eps'])
end
plot(B1g_show(jdebt),yl(1),'kx','linewidth',1.5)
plot(B1g_show(jdebt)*ones(2,1),[yl(1),PresentValueCurve(jdebt)],'k:','linewidth',1)
plot([B1g_show(jdebt2)+.001 B1g_show(jdebt)-.003],ones(2,1) * (PresentValueCurve(jdebt)+yl(1))/2, '-', 'color',[.1 .1 .1])
plot(B1g_show(jdebt)-.003,(PresentValueCurve(jdebt)+yl(1))/2, '>', 'color',[.1 .1 .1])
if savediag
	eval(['print -dpng ../../Graphs/2periods/diagram_nodefault2.png -r300'])
	eval(['print -depsc ../../Graphs/2periods/diagram_nodefault2.eps'])
end
plot(B1g_show, FundingCurve,'-.','linewidth',1,'color',PLOT.red)
plot(B1g_show(jdebt2),FundingCurve(jdebt2)-0.015,'^','color',[.1 .1 .1])
plot(B1g_show(jdebt2)*ones(2,1),[FundingCurve2(jdebt2)+0.015,FundingCurve(jdebt2)-0.015],'color',[.1 .1 .1])

if savediag
	eval(['print -dpng ../../Graphs/2periods/diagram_nodefault3.png -r300'])
	eval(['print -depsc ../../Graphs/2periods/diagram_nodefault3.eps'])
end

% Find the policy that implements Ybar in period 1
if max(N1_show) < 1 || min(N1_show) > 1
	jpot = [];
else jpot = find(N1_show >= 1, 1);
end

figure('position',[400   300   900   600])
subplot(221), hold on, grid on, ylabel('% of steady-state GDP')
plot(B1g_show,	100 * C1s_show,	'linewidth',1.5)
plot(B1g_show,	100 * C1b_show,	'linewidth',1.5)
plot(B1g_show,	100 * N1_show,	'linewidth',1.5)
l = legend('$C_1^s$','$C_1^b$','$Y$','location','southeast');
set(l,'interpreter','latex');
xlabel('$B_1^g/\bar{Y}$','Interpreter','Latex'), ylabel('% of steady-state GDP')
xlim([min(B1g_show), max(B1g_show)])
temp = get(gca,'YLim');
title('Equilibrium, t = 1')

subplot(222), hold on, grid on, ylabel('% of steady-state GDP')
plot(B1g_show,	100 * C2s_show,	'linewidth',1.5)
plot(B1g_show,	100 * C2b_show,	'linewidth',1.5)
plot(B1g_show,	100 * N2_show,	'linewidth',1.5)
l = legend('$C_2^s$','$C_2^b$','$Y$');
set(l,'interpreter','latex');
xlabel('$B_1^g/\bar{Y}$','Interpreter','Latex'), ylabel('% of steady-state GDP')
xlim([min(B1g_show), max(B1g_show)])
ylim(temp)
title('Equilibrium, t = 2')

subplot(223), hold on, grid on
plot(B1g_show, Us,	'linewidth',1.5)
plot(B1g_show, Ub,	'linewidth',1.5)
legend('Savers','Borrowers','location','southeast')
plot(B1g_show(jpot), Us(jpot), 'k*','linewidth',1.5)
plot(B1g_show(jdebt), Us(jdebt), 'ks','linewidth',1.5)
plot(B1g_show(jneutral), Us(jneutral),'ko','linewidth',1.5)
plot(B1g_show(jpot), Ub(jpot), 'k*','linewidth',1.5)
plot(B1g_show(jdebt), Ub(jdebt), 'ks','linewidth',1.5)
plot(B1g_show(jneutral), Ub(jneutral),'ko','linewidth',1.5)
xlabel('$B_1^g/\bar{Y}$','Interpreter','Latex')
xlim([min(B1g_show), max(B1g_show)])
title('Welfare')

subplot(224), hold on, grid on
plot(B1g_show, C1r,	'linewidth',1.5)
plot(B1g_show, Sr,	'linewidth',1.5,'color',PLOT.dgreen)
l = legend('$C/Y^d$', '$S/Y^d$');
set(l,'interpreter','latex');
xlabel('$B_1^g/\bar{Y}$','Interpreter','Latex')
xlim([min(B1g_show), max(B1g_show)])
title('Savers, t = 1')

if savegraph
	eval(['print -dpng ../../Graphs/2periods/',plotname,'.png -r300'])
	eval(['print -depsc ../../Graphs/2periods/',plotname,'.eps'])
end
