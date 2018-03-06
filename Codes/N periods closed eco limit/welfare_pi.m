%% This file produces pictures detailing parts of the savers welfare
% Francisco Roldán, NYU. October 2016
clear all
close all
clc
global pi_jump delev_start
%%% Remember to set run_outside = 1 in comparetiming.m before running %%%
change_gamma = 1;
show_horizon = 64;
% First set the period at which π becomes positive
pi_jump = 5;
mid_graphs = [pi_jump, pi_jump + 8, 20]; save_mid_graphs = 0;
graph_end = 1;
makedefaultgraph = 1;

g_A = 10;
pi_vec = [0 1];
delev_start_vec = pi_jump:20; % 20
Vs_compare = NaN(length(pi_vec),3,length(delev_start_vec));
Vb_compare = Vs_compare;
Cs_equiv_comp = Vs_compare;
Cb_equiv_comp = Vs_compare;
Y_lost_beta_delev_comp = NaN(length(pi_vec),1,length(delev_start_vec));
Y_lost_beta_risk_comp = NaN(length(pi_vec),1,length(delev_start_vec));
Bgvec_comp = NaN(length(delev_start_vec),show_horizon);
for jpi = 1:length(pi_vec)
	pi_crisis = pi_vec(jpi);
	for jdelev = 1:length(delev_start_vec)
		delev_start = delev_start_vec(jdelev);

		main_code_N
		if variable_speed
			pathname = 'variable_speed/';
		else
			pathname = 'parallel_paths/';
		end
		% Save results
		%%%% With EZ, welfare is consumption
		Vs_compare(jpi,:,jdelev) = 100*(Vs_mat(:) - Vs_mat(1))/Vs_mat(1); % as difference from no default benchmark
		Vb_compare(jpi,:,jdelev) = 100*(Vb_mat(:) - Vb_mat(1))/Vb_mat(1);
		% Y_lost_comp(jdelev) = Y_lost;
		Y_lost_beta_delev_comp(jpi,jdelev) = 100 * Y_lost_beta_delev; % To have NPV of losses in percent of GDP, not quarterly GDP
		Y_lost_beta_risk_comp(jpi,jdelev) = 100 * Y_lost_beta_risk;

		if pi_crisis
			Bgvec_comp(jdelev,:) = Bg_mat(3,1:show_horizon);
		end

		if ~isempty(find(mid_graphs == delev_start,1)) && save_mid_graphs
			plot_results
			if pi_crisis == 1
				pin = '_crisis';
			else pin = '_normal';
			end
			figname = ['pi',pin,'delev_start',num2str(delev_start)]
			eval(['print -dpng ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/',figname,'.png -r300'])
			eval(['print -depsc ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/',figname,'.eps'])
			% figure(1)
			% eval(['print -dpng ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/rets',figname,'.png -r300'])
		end
	end
end

delev_start_vec_plot = delev_start_vec - pi_jump;
% set(0,'defaulttextinterpreter','latex')
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(groot, 'defaultLegendInterpreter','latex');
lw = 1.5;


figure, title('Output losses (no default path)'), ylabel('% of steady-state GDP'), xlabel('Deleveraging delay')
hold on, grid on
p1 = plot(delev_start_vec_plot,Y_lost_beta_delev_comp(1,:),'linewidth',lw,'color',PLOT.blue);
% plot(delev_start_vec_plot,Y_lost_beta_risk_comp(1,:),'--','linewidth',lw,'color',PLOT.blue)
p2 = plot(delev_start_vec_plot,Y_lost_beta_delev_comp(2,:),'linewidth',lw,'color',PLOT.red);
% plot(delev_start_vec_plot,Y_lost_beta_risk_comp(2,:),'--','linewidth',lw,'color',PLOT.red)
l = legend([p1 p2],['$\pi = ','\pi^{normal}','$'],['$\pi = ','\pi^{crisis}','$']);
set(l,'interpreter','latex','location','southeast')
if graph_end
	eval(['print -dpng ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/outputloss_pi.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/outputloss_pi.eps'])
end

figure('position', [100 300 1250 550])
subplot(131), grid on, hold on, xlabel('Deleveraging delay'), title('Consumption equivalents - Savers'), xlim([0 14]), ylabel('% of permanent consumption without risk')
p = Vs_compare(2,2,:);	% no delev, pi = pi_crisis
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.blue)
p = Vs_compare(2,3,:);  % delev, pi = pi_crisis
p1 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.red);
p = Vs_compare(2,1,:);	% no risk, pi = pi_crisis
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.yellow)
legend('no deleveraging','deleveraging', 'no risk')
subplot(132), grid on, hold on, xlabel('Deleveraging delay'), title('Consumption equivalents - Borrowers'), xlim([0 14]), ylabel('% of permanent consumption without risk')
p = Vb_compare(2,2,:);	% no delev, pi = pi_crisis
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.blue)
p = Vb_compare(2,3,:);  % delev, pi = pi_crisis
p1 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.red);
p = Vb_compare(2,1,:);	% no risk, pi = pi_crisis
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.yellow)
legend('no deleveraging','deleveraging', 'no risk','location','northeast')

subplot(133), hold on, grid on, ylabel('% of steady-state GDP'), xlabel('Deleveraging delay'), title('Output losses (no default path)'), xlim([0 14])
plot(delev_start_vec_plot,Y_lost_beta_delev_comp(2,:),'linewidth',lw,'color',PLOT.blue)
% plot(delev_start_vec_plot,Y_lost_beta_risk_comp(2,:),'--','linewidth',lw,'color',PLOT.blue)
l = legend(['$\pi = ','\pi^{crisis}','$']);
set(l,'interpreter','latex')
if graph_end
	eval(['print -dpng ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/timing_comparison.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/timing_comparison.eps'])
end

figure('position',[700 500 900 400])
subplot(121), hold on, grid on, xlabel('Deleveraging delay'), xlim([0 14])
p = Vs_compare(1,2,:);	% gamma = 0.1, no delev
% plot(delev_start_vec_plot,p(:),'--','linewidth',lw,'color',PLOT.blue)
p = Vs_compare(1,3,:);	% gamma = 0.1, delev
p1 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.blue);
p = Vs_compare(1,1,:);	% gamma = 0.1, no risk
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.lgray)
p = Vs_compare(2,2,:);	% gamma = 1, no delev
% plot(delev_start_vec_plot,p(:),'--','linewidth',lw,'color',PLOT.red)
p = Vs_compare(2,3,:);	% gamma = 1, delev
p2 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.red);
title('Consumption equivalents - Savers'), ylabel('% of permanent consumption without risk')
l = legend([p1 p2],['$\pi = ','\pi^{normal}','$'],['$\pi = ','\pi^{crisis}','$']);
set(l, 'interpreter','latex','location','east')

subplot(122), hold on, grid on, xlabel('Deleveraging delay'), xlim([0 14])
p = Vb_compare(1,2,:); % V(gamma, risk/delev, delay)
% plot(delev_start_vec_plot,p(:),'--','linewidth',lw,'color',PLOT.blue)
p = Vb_compare(1,3,:);
p1 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.blue);
p = Vb_compare(1,1,:);
plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.lgray)
p = Vb_compare(2,2,:); % V(gamma, risk/delev, delay)
% plot(delev_start_vec_plot,p(:),'--','linewidth',lw,'color',PLOT.red)
p = Vb_compare(2,3,:);
p2 = plot(delev_start_vec_plot,p(:),'linewidth',lw,'color',PLOT.red);
title('Consumption equivalents - Borrowers'), ylabel('% of permanent consumption without risk')
l = legend([p1 p2],['$\pi = ','\pi^{normal}','$'],['$\pi = ','\pi^{crisis}','$']);
set(l, 'interpreter','latex','location','east')
if graph_end
	eval(['print -dpng ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/ConsEquiv_pi.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/ConsEquiv_pi.eps'])
end

Pbvec = linspace(.5,1,1000);
opt_d = NaN(length(pi_vec),length(Pbvec));
for jpi = 1:length(pi_vec)
	for jp = 1:length(Pbvec)
		Pb = Pbvec(jp); Ps = 1-Pb;
		welf(:) = Pb * Vb_compare(jpi,3,:) + (Ps) * Vs_compare(jpi,3,:);
		[w,opt_d(jpi,jp)] = max(welf);
	end
end
figure, hold on, grid on, title('Optimal delay'), xlabel('Relative Pareto weight of borrowers')
plot(Pbvec, opt_d(1,:),'linewidth',1.5,'color',PLOT.blue)
plot(Pbvec, opt_d(2,:),'linewidth',1.5,'color',PLOT.red)
l = legend(['$\pi = ','\pi^{normal}','$'],['$\pi = ','\pi^{crisis}','$'],'location','southeast');
set(l,'interpreter','latex')
if graph_end
	eval(['print -dpng ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/optimaldelay.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/optimaldelay.eps'])
end


figure('position',[700 479 668 421])
subplot(121), hold on, grid on,	title('Household debt')
plot(1:show_horizon,100*chi*Bhvec_mat(1,1:show_horizon)/ YearGDP,'linewidth',lw)
xlabel('Quarters'), ylabel('% of steady-state GDP')
subplot(122), hold on, grid on, title('Government debt')
plot(1:show_horizon,100*Bgvec_comp(1,1:show_horizon)/ YearGDP,'linewidth',lw)
plot(1:show_horizon,100*Bgvec_comp(8,1:show_horizon)/ YearGDP,'linewidth',lw)
plot(1:show_horizon,100*Bgvec_comp(end,1:show_horizon)/ YearGDP,'linewidth',lw)
xlabel('Quarters'), ylabel('% of steady-state GDP')
legend('Immediate','Half','Late')
if graph_end
	eval(['print -dpng ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/delevpath.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/delevpath.eps'])
end

Bgpath = Bgvec_comp(1,1:show_horizon);
Bgdef1 = Bgpath; Bgdef1(18:end) = (1-haircut) * Bgpath(18:end); Bgdef1(1:14) = NaN;
Bgdef2 = Bgpath; Bgdef2(26:end) = (1-haircut) * Bgpath(26:end); Bgdef2(1:22) = NaN;
figure('position',[680 558 560 420]), hold on, grid on, title('Government debt'), ylabel('% of steady-state GDP'), xlabel('Quarters')
plot(1:show_horizon,100*Bgpath / YearGDP,'linewidth',lw,'color',PLOT.blue)
plot(1:show_horizon,100*Bgdef1 / YearGDP,'--','linewidth',lw,'color',PLOT.red)
plot(1:show_horizon,100*Bgdef2 / YearGDP,':','linewidth',lw,'color',PLOT.yellow)
ylim([30, 120])
if makedefaultgraph
	eval(['print -dpng ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/debtpath.png -r300'])
	eval(['print -depsc ../../Graphs/Nperiods/Constant_tax_rate/',pathname,'pi/debtpath.eps'])
end
