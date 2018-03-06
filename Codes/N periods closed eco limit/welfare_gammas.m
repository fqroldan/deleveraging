%% This file produces pictures detailing parts of the savers welfare
% Francisco Roldán, NYU. 2016
clear all
close all
clc
global pi_jump delev_start welfare_gammas_running highdebt
welfare_gammas_running = 1;
change_chi = 1; % benchmark = 1

OpenEco = 0;
highdebt = 0;

quick_graphs = 0;

if change_chi ~= 0
	if change_chi == 1
		chi_up = 1/2;
		chiup = 'half';
		R_borr = 12;
	elseif change_chi == 2
		chi_up = 1/3;
		chiup = 'third';
		R_borr = 15;
	elseif change_chi == 3
		chi_up = 1/2;
		chiup = 'half';
		R_borr = 15;
	else
		disp('Error: change_chi must be 0, 1, 2, or 3')
		stop
	end
beta_borr_up = (1 + R_borr/100)^(-1/4);	% deduce β_borr from annualized interest rate
end

%%% Remember to set run_outside = 1 in main_code_N.m before running %%%
change_gamma = 1;
show_horizon = 64;	%64
% First set the period at which π becomes positive
pi_jump = 5;
mid_graphs = [pi_jump, pi_jump + 8, 20]; save_mid_graphs = 1;
graph_end = 0;
makedefaultgraph = 0;

pi_normal = 0;

if pi_normal
	pi_crisis = 0;
end
if highdebt
	initBg = 1.2;
end

if quick_graphs
	gamma_vec = [10];
	delev_start_vec = mid_graphs;
	save_mid_graphs = 1;
else
	gamma_vec = [1 10 25];
	delev_start_vec = pi_jump:1:20; % 20
end
Vs_compare = NaN(length(gamma_vec),3,length(delev_start_vec));
Vb_compare = Vs_compare;
Cs_equiv_comp = Vs_compare;
Cb_equiv_comp = Vs_compare;
Y_lost_beta_delev_comp = NaN(length(gamma_vec),1,length(delev_start_vec));
Y_lost_beta_risk_comp = NaN(length(gamma_vec),1,length(delev_start_vec));
Bgvec_comp = NaN(length(delev_start_vec),show_horizon);
Bhvec_comp = Bgvec_comp;
q_comp = NaN(show_horizon+1,length(gamma_vec),length(delev_start_vec));
q0_comp = NaN(show_horizon+1,length(delev_start_vec));

for jg = 1:length(gamma_vec)
	g_A = gamma_vec(jg);
	for jdelev = 1:length(delev_start_vec)
		delev_start = delev_start_vec(jdelev);

		main_code_N
		if variable_speed
			pathname = 'variable_speed/';
		else
			pathname = 'parallel_paths/';
		end
		if pi_normal
			pathname = [pathname, 'pi/'];
		end
		if change_chi
			pathname = [pathname, 'chi', chiup, 'Rb', num2str(R_borr), '/'];
		end
		if OpenEco
			pathname = ['OpenEco/', pathname];
		end
		if highdebt
			pathname = [pathname, 'highBg/'];
		end
		% Save results
		%%%% With EZ, welfare is consumption
		Vs_compare(jg,:,jdelev) = 100*(Vs_mat(:) - Vs_mat(1))/Vs_mat(1); % as difference from no default benchmark
		Vb_compare(jg,:,jdelev) = 100*(Vb_mat(:) - Vb_mat(1))/Vb_mat(1);

		Cs_equiv_comp(jg,:,jdelev) = 100 * (C_equiv_s(:) - C_equiv_s(1))/C_equiv_s(1);
		Cb_equiv_comp(jg,:,jdelev) = 100 * (C_equiv_b(:) - C_equiv_b(1))/C_equiv_b(1);
		% Y_lost_comp(jdelev) = Y_lost;
		Y_lost_beta_delev_comp(jg,jdelev) = 100 * Y_lost_beta_delev; % To have NPV of losses in percent of GDP, not quarterly GDP
		Y_lost_beta_risk_comp(jg,jdelev) = 100 * Y_lost_beta_risk;

		q_comp(:,jg,jdelev) = q_mat(3,1:show_horizon+1);
		q0_comp(:,jg,jdelev) = q0_mat(3,1:show_horizon+1);

		if g_A == 10
			Bgvec_comp(jdelev,:) = Bg_mat(3,1:show_horizon);
			Bhvec_comp(jdelev,:) = Bh_mat(3,1:show_horizon);
		end

		if ~isempty(find(mid_graphs == delev_start,1)) && save_mid_graphs
			plot_results
			% figure(2)
			figname = ['gamma',num2str(g_A),'delev_start',num2str(delev_start)]
			% eval(['print -dpng ../../Graphs/Nperiods/',pathname,figname,'.png -r300'])
			% eval(['print -depsc ../../Graphs/Nperiods/',pathname,figname,'.eps'])
			% figure(3)
			% eval(['print -dpng ../../Graphs/Nperiods/',pathname,'debts',figname,'.png -r300'])
			% eval(['print -depsc ../../Graphs/Nperiods/',pathname,'debts',figname,'.eps'])
			figure(4)
			eval(['print -dpng ../../Graphs/Nperiods/',pathname,'sum',figname,'.png -r300'])
			eval(['print -depsc ../../Graphs/Nperiods/',pathname,'sum',figname,'.eps'])
		end
	end
end
close all

delev_start_vec_plot = delev_start_vec - pi_jump;
% set(0,'defaulttextinterpreter','latex')
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(groot, 'defaultLegendInterpreter','latex');
lw = 1.5;


figure(1), title('Capitalized output losses (no default path)'), ylabel('% of steady-state GDP'), xlabel('Deleveraging delay')
hold on, grid on
p1 = plot(delev_start_vec_plot,Y_lost_beta_delev_comp(1,:),'linewidth',lw,'color',PLOT.blue);
% plot(delev_start_vec_plot,Y_lost_beta_risk_comp(1,:),'--','linewidth',lw,'color',PLOT.blue)
p2 = plot(delev_start_vec_plot,Y_lost_beta_delev_comp(2,:),'linewidth',lw,'color',PLOT.red);
% plot(delev_start_vec_plot,Y_lost_beta_risk_comp(2,:),'--','linewidth',lw,'color',PLOT.red)
p3 = plot(delev_start_vec_plot,Y_lost_beta_delev_comp(3,:),'linewidth',lw,'color',PLOT.dgreen);
% plot(delev_start_vec_plot,Y_lost_beta_risk_comp(3,:),'--','linewidth',lw,'color',PLOT.dgreen)
% p4 = plot(delev_start_vec_plot,Y_lost_beta_delev_comp(4,:),'linewidth',lw,'color',PLOT.dgreen);
% plot(delev_start_vec_plot,Y_lost_beta_risk_comp(4,:),'--','linewidth',lw,'color',PLOT.dgreen)
l = legend([p1 p2 p3],['$\gamma = ',num2str(gamma_vec(1)),'$'],['$\gamma = ',num2str(gamma_vec(2)),'$'],['$\gamma = ',num2str(gamma_vec(3)),'$']);
set(l,'interpreter','latex','location','southeast')
if graph_end
	eval(['print -dpng ../../Graphs/Nperiods/',pathname,'outputloss_gammas.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/',pathname,'outputloss_gammas.eps'])
end

figure(2)
set(2,'position', [100 300 1250 450])
subplot(131), grid on, hold on, xlabel('Deleveraging delay'), title('Consumption equivalents - Savers'), xlim([0 14]), ylabel('% of permanent consumption without risk')
% p = Cs_equiv_comp(2,2,:);	% gamma = 10, no delev
% plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.blue)
p = Cs_equiv_comp(2,3,:);  % gamma = 10, delev
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.red);
p = Cs_equiv_comp(2,1,:);	% gamma = 10, no risk
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.yellow)
% legend('no deleveraging','deleveraging', 'no risk')
legend('deleveraging', 'no risk')
subplot(132), grid on, hold on, xlabel('Deleveraging delay'), title('Consumption equivalents - Borrowers'), xlim([0 14]), ylabel('% of permanent consumption without risk')
% p = Cb_equiv_comp(2,2,:);	% gamma = 10, no delev
% plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.blue)
p = Cb_equiv_comp(2,3,:);  % gamma = 10, delev
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.red);
p = Cb_equiv_comp(2,1,:);	% gamma = 10, no risk
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.yellow)
% legend('no deleveraging','deleveraging', 'no risk','location','northeast')
legend('deleveraging', 'no risk','location','northeast')

subplot(133), hold on, grid on, ylabel('% of steady-state GDP'), xlabel('Deleveraging delay'), title('Capitalized output losses (no default path)'), xlim([0 14])
plot(delev_start_vec_plot,Y_lost_beta_delev_comp(2,:),'linewidth',lw,'color',PLOT.blue)
% plot(delev_start_vec_plot,Y_lost_beta_risk_comp(2,:),'--','linewidth',lw,'color',PLOT.blue)
l = legend(['$\gamma = ',num2str(gamma_vec(2)),'$']);
set(l,'interpreter','latex')
if graph_end
	eval(['print -dpng ../../Graphs/Nperiods/',pathname,'timing_comparison.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/',pathname,'timing_comparison.eps'])
end

figure(3)
set(3,'position',[584 232 950 450])
subplot(121), hold on, grid on, xlabel('Deleveraging delay'), xlim([0 14])
p = Cs_equiv_comp(1,2,:);	% gamma = 1, no delev
% plot(delev_start_vec_plot,p(:),'--','linewidth',lw,'color',PLOT.blue)
p = Cs_equiv_comp(1,3,:);	% gamma = 1, delev
p1 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.blue);
p = Cs_equiv_comp(1,1,:);	% gamma = 0.1, no risk
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.lgray)
p = Cs_equiv_comp(2,2,:);	% gamma = 10, no delev
% plot(delev_start_vec_plot,p(:),'--','linewidth',lw,'color',PLOT.red)
p = Cs_equiv_comp(2,3,:);	% gamma = 10, delev
p2 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.red);
p = Cs_equiv_comp(3,2,:);	% gamma = 100, no delev
% plot(delev_start_vec_plot,p(:),'--','linewidth',lw,'color',PLOT.dgreen)
p = Cs_equiv_comp(3,3,:);	% gamma = 100, delev
p3 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.dgreen);
title('Consumption equivalents - Savers'), ylabel('% of permanent consumption without risk')
l = legend([p1 p2 p3],['$\gamma = ',num2str(gamma_vec(1)),'$'],['$\gamma = ',num2str(gamma_vec(2)),'$'],['$\gamma = ',num2str(gamma_vec(3)),'$']);
set(l, 'interpreter','latex','location','east')

subplot(122), hold on, grid on, xlabel('Deleveraging delay'), xlim([0 14])
p = Cb_equiv_comp(1,2,:); % V(gamma, risk/delev, delay)
% plot(delev_start_vec_plot,p(:),'--','linewidth',lw,'color',PLOT.blue)
p = Cb_equiv_comp(1,3,:);
p1 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.blue);
p = Cb_equiv_comp(1,1,:);
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.lgray)
p = Cb_equiv_comp(2,2,:); % V(gamma, risk/delev, delay)
% plot(delev_start_vec_plot,p(:),'--','linewidth',lw,'color',PLOT.red)
p = Cb_equiv_comp(2,3,:);
p2 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.red);
p = Cb_equiv_comp(3,2,:);
% plot(delev_start_vec_plot,p(:),'--','linewidth',lw,'color',PLOT.dgreen)
p = Cb_equiv_comp(3,3,:);
p3 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.dgreen);
title('Consumption equivalents - Borrowers'), ylabel('% of permanent consumption without risk')
l = legend([p1 p2 p3],['$\gamma = ',num2str(gamma_vec(1)),'$'],['$\gamma = ',num2str(gamma_vec(2)),'$'],['$\gamma = ',num2str(gamma_vec(3)),'$']);
set(l, 'interpreter','latex','location','east')
if graph_end
	eval(['print -dpng ../../Graphs/Nperiods/',pathname,'ConsEquiv_gammas.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/',pathname,'ConsEquiv_gammas.eps'])
end

Pbvec = linspace(.5,1,1000);
opt_d = NaN(length(gamma_vec),length(Pbvec));
for jg = 1:length(gamma_vec)
	for jp = 1:length(Pbvec)
		Pb = Pbvec(jp); Ps = 1-Pb;
		welf(:) = Pb * Vb_compare(jg,3,:) + (Ps) * Vs_compare(jg,3,:);
		[w,opt_d(jg,jp)] = max(welf);
	end
end
figure(4), hold on, grid on, title('Optimal delay'), xlabel('Relative Pareto weight of borrowers')
plot(Pbvec, opt_d(1,:),'linewidth',1.5,'color',PLOT.blue)
plot(Pbvec, opt_d(2,:),'linewidth',1.5,'color',PLOT.red)
plot(Pbvec, opt_d(3,:),'linewidth',1.5,'color',PLOT.dgreen)
l = legend(['$\gamma = ',num2str(gamma_vec(1)),'$'],['$\gamma = ',num2str(gamma_vec(2)),'$'],['$\gamma = ',num2str(gamma_vec(3)),'$'],'location','southeast');
set(l,'interpreter','latex')
if graph_end
	eval(['print -dpng ../../Graphs/Nperiods/',pathname,'optimaldelay.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/',pathname,'optimaldelay.eps'])
end


figure(5)
subplot(121), hold on, grid on,	title('Household debt')
plot(1:show_horizon,100 * chi*Bhvec_comp(1,1:show_horizon)/ YearGDP,'linewidth',lw)
xlabel('Quarters'), ylabel('% of steady-state GDP')
subplot(122), hold on, grid on, title('Government debt')
plot(1:show_horizon,100 * Bgvec_comp(1,1:show_horizon)/ YearGDP,'linewidth',lw)
plot(1:show_horizon,100 * Bgvec_comp(8,1:show_horizon)/ YearGDP,'linewidth',lw)
plot(1:show_horizon,100 * Bgvec_comp(end,1:show_horizon)/ YearGDP,'linewidth',lw)
xlabel('Quarters'), ylabel('% of steady-state GDP')
legend('Early','Half','Late')
if graph_end
	eval(['print -dpng ../../Graphs/Nperiods/',pathname,'delevpath.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/',pathname,'delevpath.eps'])
end

Bgpath = Bgvec_comp(1,1:show_horizon);
Bgdef1 = Bgpath; Bgdef1(18:end) = (1-haircut) * Bgpath(18:end); Bgdef1(1:14) = NaN;
Bgdef2 = Bgpath; Bgdef2(26:end) = (1-haircut) * Bgpath(26:end); Bgdef2(1:22) = NaN;
figure(6), hold on, grid on, title('Government debt'), ylabel('% of steady-state GDP'), xlabel('Quarters')
plot(1:show_horizon,100*Bgpath / YearGDP,'linewidth',lw,'color',PLOT.blue)
plot(1:show_horizon,100*Bgdef1 / YearGDP,'--','linewidth',lw,'color',PLOT.red)
plot(1:show_horizon,100*Bgdef2 / YearGDP,':','linewidth',lw,'color',PLOT.yellow)
ylim([30, 120])
if makedefaultgraph
	eval(['print -dpng ../../Graphs/Nperiods/',pathname,'debtpath.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/',pathname,'debtpath.eps'])
end

figure(7), hold on, grid on, title('Price of govt debt')
l = [];
p1 = plot(0:show_horizon, q_comp(:,1,1),'linewidth',1.5,'color',PLOT.c{1});
plot(0:show_horizon,q0_comp(:,1,1),':','linewidth',1.5,'color',PLOT.c{1});
p2 = plot(0:show_horizon, q_comp(:,2,1),'linewidth',1.5,'color',PLOT.c{2});
plot(0:show_horizon,q0_comp(:,2,1),':','linewidth',1.5,'color',PLOT.c{2});
p3 = plot(0:show_horizon, q_comp(:,3,1),'linewidth',1.5,'color',PLOT.c{3});
plot(0:show_horizon,q0_comp(:,3,1),':','linewidth',1.5,'color',PLOT.c{3});

l = legend([p1 p2 p3], ['$\gamma = ',num2str(gamma_vec(1)),'$'],['$\gamma = ',num2str(gamma_vec(2)),'$'],['$\gamma = ',num2str(gamma_vec(3)),'$']);
xlabel('Quarters'), set(l,'interpreter','latex')
if makedefaultgraph
	eval(['print -dpng ../../Graphs/Nperiods/',pathname,'debtprice.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/',pathname,'debtprice.eps'])
end
clear welfare_gammas_running
