lw = 1.5;
le = [];
start_default_path = 30;
start_postdefault_path = min(start_default_path + floor(1/theta), show_horizon+1);

Cs1_show = 100*(Cs_mat(1,1:show_horizon+1)-Cs_SS)/Cs_SS;
Cs2_show = 100*(Cs_mat(2,1:show_horizon+1)-Cs_SS)/Cs_SS;
Cs3_show = 100*(Cs_mat(3,1:show_horizon+1)-Cs_SS)/Cs_SS;
Cs4_show = 100*(Csdi_mat(3,1:show_horizon+1)-Cs_SS)/Cs_SS;
Cs5_show = 100*(Csd_mat(3,1:show_horizon+1)-Cs_SS)/Cs_SS;
% Cs1_show(1) = Cs3_show(1); Cs2_show(1) = Cs3_show(1);
Cs4_show(1:start_default_path) = Cs3_show(1:start_default_path);
Cs4_show(1:start_default_path-4) = NaN;
Cs5_show(1:start_postdefault_path) = Cs4_show(1:start_postdefault_path);
Cs5_show(1:start_postdefault_path-4) = NaN;
% Cs4_show(start_postdefault_path:end) = Cs5_show(start_postdefault_path:end);


Cb1_show = 100*(Cb_mat(1,1:show_horizon+1)-Cb_SS)/Cb_SS;
Cb2_show = 100*(Cb_mat(2,1:show_horizon+1)-Cb_SS)/Cb_SS;
Cb3_show = 100*(Cb_mat(3,1:show_horizon+1)-Cb_SS)/Cb_SS;
Cb4_show = 100*(Cbdi_mat(3,1:show_horizon+1)-Cb_SS)/Cb_SS;
Cb5_show = 100*(Cbd_mat(3,1:show_horizon+1)-Cb_SS)/Cb_SS;
% Cb1_show(1) = Cb3_show(1); Cb2_show(1) = Cb3_show(1);
Cb4_show(1:start_default_path) = Cb3_show(1:start_default_path);
Cb4_show(1:start_default_path-4) = NaN;
Cb5_show(1:start_postdefault_path) = Cb4_show(1:start_postdefault_path);
Cb5_show(1:start_postdefault_path-4) = NaN;
% Cb4_show(start_postdefault_path:end) = Cb5_show(start_postdefault_path:end);
E2_show = E2_mat(3,1:show_horizon+1);
CbE_show = Cb3_show;
CbE_show(E2_show < 0) = NaN;
CbE_show(1:pi_jump+1) = NaN;

Y1_show = 100*(Y_mat(1,1:show_horizon+1)-Y_SS)/Y_SS;
Y2_show = 100*(Y_mat(2,1:show_horizon+1)-Y_SS)/Y_SS;
Y3_show = 100*(Y_mat(3,1:show_horizon+1)-Y_SS)/Y_SS;
Y4_show = 100*(Ydi_mat(3,1:show_horizon+1)-Y_SS)/Y_SS;
Y5_show = 100*(Yd_mat(3,1:show_horizon+1)-Y_SS)/Y_SS;
% Y1_show(1) = Y3_show(1); Y2_show(1) = Y3_show(1);
Y4_show(1:start_default_path) = Y3_show(1:start_default_path);
Y4_show(1:start_default_path-4) = NaN;
Y5_show(1:start_postdefault_path) = Y4_show(1:start_postdefault_path);
Y5_show(1:start_postdefault_path-4) = NaN;
% Y4_show(start_postdefault_path:end) = Y5_show(start_postdefault_path:end);

pi1_show = (1+pidef_mat(1,1:show_horizon+1)).^(1+3*(1-annual))-1;
pi2_show = (1+pidef_mat(2,1:show_horizon+1)).^(1+3*(1-annual))-1;
pi3_show = (1+pidef_mat(3,1:show_horizon+1)).^(1+3*(1-annual))-1;
% pi1_show(1) = pi3_show(1); pi2_show(1) = pi3_show(1);

q1_show = q_mat(1,1:show_horizon+1);
q2_show = q_mat(2,1:show_horizon+1);
q3_show = q_mat(3,1:show_horizon+1);
% q1_show(1) = q3_show(1); q2_show(1) = q3_show(1);
q0_show = q0_mat(3,1:show_horizon+1);

T1_show = 100*(T_mat(1,1:show_horizon+1))/Y_SS;
T2_show = 100*(T_mat(2,1:show_horizon+1))/Y_SS;
T3_show = 100*(T_mat(3,1:show_horizon+1))/Y_SS;
T4_show = 100*(Tdi_mat(3,1:show_horizon+1))/Y_SS;
T5_show = 100*(Td_mat(3,1:show_horizon+1))/Y_SS;
% T1_show(1) = T3_show(1); T2_show(1) = T3_show(1);
T4_show(1:start_default_path) = T3_show(1:start_default_path);
T4_show(1:start_default_path-4) = NaN;
T5_show(1:start_postdefault_path) = T4_show(1:start_postdefault_path);
T5_show(1:start_postdefault_path-4) = NaN;
% T4_show(start_postdefault_path:end) = T5_show(start_postdefault_path:end);

Tr1_show = 100 * T_mat(1,1:show_horizon+1)./Y_mat(1,1:show_horizon+1);
Tr2_show = 100 * T_mat(2,1:show_horizon+1)./Y_mat(2,1:show_horizon+1);
Tr3_show = 100 * T_mat(3,1:show_horizon+1)./Y_mat(3,1:show_horizon+1);
Tr4_show = 100 * Tdi_mat(3,1:show_horizon+1)./Ydi_mat(3,1:show_horizon+1);
Tr5_show = 100 * Td_mat(3,1:show_horizon+1)./Yd_mat(3,1:show_horizon+1);
% Tr1_show(1) = Tr3_show(1); Tr2_show(1) = Tr3_show(1);
Tr4_show(1:start_default_path) = Tr3_show(1:start_default_path);
Tr4_show(1:start_default_path-4) = NaN;
Tr5_show(1:start_postdefault_path) = Tr4_show(1:start_postdefault_path);
Tr5_show(1:start_postdefault_path-4) = NaN;
% Tr4_show(start_postdefault_path:end) = Tr5_show(start_postdefault_path:end);

BY1_show = 100 * Bg_mat(1,1:show_horizon+1)./(4*Y_mat(1,1:show_horizon+1));
BY2_show = 100 * Bg_mat(2,1:show_horizon+1)./(4*Y_mat(2,1:show_horizon+1));
BY3_show = 100 * Bg_mat(3,1:show_horizon+1)./(4*Y_mat(3,1:show_horizon+1));
BY4_show = 100 * (1-haircut) * Bg_mat(3,1:show_horizon+1) ./ (4*Ydi_mat(3,1:show_horizon+1));
BY4_show(1:start_default_path) = BY3_show(1:start_default_path);
BY4_show(1:start_default_path-4) = NaN;

%
% figure(1), set(gcf,'position', [0, 150, 1000, 750])
% subplot(121), hold on, grid on, xlabel('Quarters'), ylabel('% of GDP'), xlim([0 5*ceil(show_horizon/5)])
% 	for j = 1:3
% 		plot(0:show_horizon-1,Bhvec_mat(j,1:show_horizon)*chi / (N(end)*(1+3*(1-annual))), 'linewidth',lw);
% 	end
% 	title('Household debt')
% subplot(122), hold on, grid on, xlabel('Quarters'), ylabel('% of GDP'), xlim([0 5*ceil(show_horizon/5)])
% 	for j = 1:3
% 		plot(0:show_horizon-1,Bgvec_mat(j,1:show_horizon)/(N(end)*(1+3*(1-annual))), 'linewidth',lw);
% 	end
% 	title('Government debt'), legend('no deleveraging', 'deleveraging', 'no risk')

if plot_returns
	figure(1)
	subplot(121), hold on, grid on, title('(Risk-neutral) Returns')
	p1 = plot(1:show_horizon, ret_0_mat(2,1:show_horizon)-1, 'linewidth',lw,'color',PLOT.blue);
	plot(1:show_horizon, ret_def_mat(2,1:show_horizon)-1, '--', 'linewidth',lw, 'color',PLOT.blue)
	p2 = plot(1:show_horizon, ret_0_mat(3,1:show_horizon)-1, 'linewidth',lw,'color',PLOT.red);
	plot(1:show_horizon, ret_def_mat(3,1:show_horizon)-1, '--', 'linewidth',lw, 'color',PLOT.red)
	p3 = plot(1:show_horizon, ret_0_mat(1,1:show_horizon)-1, 'linewidth',lw,'color',PLOT.yellow);
	plot(1:show_horizon, ret_def_mat(1,1:show_horizon)-1, '--', 'linewidth',lw, 'color',PLOT.yellow)
	legend([p1 p2 p3],'no deleveraging','deleveraging','no risk','location','east')
	xlabel('Quarters'), xlim([0 show_horizon])
	subplot(122), hold on, grid on, title('(EZ) Returns')
	p1 = plot(1:show_horizon, ret_EZ_0_mat(2,1:show_horizon)-1, 'linewidth',lw,'color',PLOT.blue);
	plot(1:show_horizon, ret_EZ_def_mat(2,1:show_horizon)-1, '--', 'linewidth',lw, 'color',PLOT.blue)
	p2 = plot(1:show_horizon, ret_EZ_0_mat(3,1:show_horizon)-1, 'linewidth',lw,'color',PLOT.red);
	plot(1:show_horizon, ret_EZ_def_mat(3,1:show_horizon)-1, '--', 'linewidth',lw, 'color',PLOT.red)
	p3 = plot(1:show_horizon, ret_EZ_0_mat(1,1:show_horizon)-1, 'linewidth',lw,'color',PLOT.yellow);
	plot(1:show_horizon, ret_EZ_def_mat(1,1:show_horizon)-1, '--', 'linewidth',lw, 'color',PLOT.yellow)
	legend([p1 p2 p3],'no deleveraging','deleveraging','no risk','location','east')
	xlabel('Quarters'), xlim([0 show_horizon])
end
disp('graphs')
figure(2), set(gcf,'position', [200 150 1000 750])
subplot(2,3,1), hold on, grid on, xlabel('Quarters'), ylabel('%'), xlim([0 5*ceil(show_horizon/5)])
	plot(0:show_horizon,Cs1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% plot(0:show_horizon,Cs2_show,'-.','linewidth',lw,'color',PLOT.blue);
	plot(0:show_horizon,Cs3_show,'linewidth',lw,'color',PLOT.red);
	plot(0:show_horizon,Cs4_show,':','color',PLOT.purple,'linewidth',lw);
	plot(0:show_horizon,Cs5_show,':','color',PLOT.blue,'linewidth',lw);
	% plot(0:show_horizon,Cs5_show,':','color',PLOT.dgreen,'linewidth',lw);
	title('Savers consumption')
subplot(2,3,2), hold on, grid on, xlabel('Quarters'), ylabel('%'), xlim([0 5*ceil(show_horizon/5)])
	plot(0:show_horizon,Cb1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% plot(0:show_horizon,Cb2_show,'-.','linewidth',lw,'color',PLOT.blue);
	plot(0:show_horizon,Cb3_show,'linewidth',lw,'color',PLOT.red);
	plot(0:show_horizon,Cb4_show,':','color',PLOT.purple,'linewidth',lw);
	plot(0:show_horizon,Cb5_show,':','color',PLOT.blue,'linewidth',lw);
	plot(0:show_horizon,CbE_show,'*','linewidth',lw,'color',PLOT.red);
	% set(gca, 'Ylim', [-30 5])
	% yl = get(gca,'Ylim'); E2_show = E2_show - (E2_show - min(yl) - .1) .* (E2_show == min(E2_show));
	% area(0:show_horizon,E2_show * (max(yl)-1e-1),min(yl)+1e-1,'facecolor',[.95 .85 .85],'linestyle','none');
	% set(gca,'Ylim',yl);
	% plot(0:show_horizon,Cb1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% plot(0:show_horizon,Cb2_show,'-.','linewidth',lw,'color',PLOT.blue);
	% plot(0:show_horizon,Cb3_show,'linewidth',lw,'color',PLOT.red);
	% plot(0:show_horizon,Cb4_show,':','color',PLOT.purple,'linewidth',lw);
	title('Borrowers consumption')
subplot(2,3,3), hold on, grid on, xlabel('Quarters'), ylabel('%'), xlim([0 5*ceil(show_horizon/5)])
	p1 = plot(0:show_horizon,Y1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% p2 = plot(0:show_horizon,Y2_show,'-.','linewidth',lw,'color',PLOT.blue);
	p3 = plot(0:show_horizon,Y3_show,'linewidth',lw,'color',PLOT.red);
	p4 = plot(0:show_horizon,Y4_show,':','color',PLOT.purple,'linewidth',lw);
	plot(0:show_horizon,Y5_show,':','color',PLOT.blue,'linewidth',lw);
	title('Output')
	% legend([p1 p2 p3], 'no sovereign risk', 'no deleveraging', 'deleveraging', 'location', 'southeast')
	legend([p1 p3], 'no sovereign risk', 'deleveraging', 'location', 'southeast')
	% set(gca, 'Ylim', [-14 2])
subplot(2,3,4), hold on, grid on, xlabel('Quarters'), xlim([0 5*ceil(show_horizon/5)])
	plot(0:show_horizon,pi1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% plot(0:show_horizon,pi2_show,'-.','linewidth',lw,'color',PLOT.blue);
	plot(0:show_horizon,pi3_show,'linewidth',lw,'color',PLOT.red);
	title('Default probability (yearly)')
subplot(2,3,5), hold on, grid on, xlabel('Quarters'), xlim([0 5*ceil(show_horizon/5)])
	plot(0:show_horizon,q1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% plot(0:show_horizon,q2_show,'-.','linewidth',lw,'color',PLOT.blue);
	plot(0:show_horizon,q3_show,'linewidth',lw,'color',PLOT.red);
	plot(0:show_horizon,q0_show,'--', 'linewidth',lw,'color',PLOT.dgreen);
	title('Price of government debt')
subplot(2,3,6), hold on, grid on, xlabel('Quarters'), ylabel('% of steady-state GDP'), xlim([0 5*ceil(show_horizon/5)])
	plot(0:show_horizon,T1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% plot(0:show_horizon,T2_show,'-.','linewidth',lw,'color',PLOT.blue);
	plot(0:show_horizon,T3_show,'linewidth',lw,'color',PLOT.red);
	plot(0:show_horizon,T4_show,':','color',PLOT.purple,'linewidth',lw);
	title('Tax Collections')

figure(4), set(gcf,'position', [200 150 1000 750])
subplot(2,2,1), hold on, grid on, xlabel('Quarters'), ylabel('%'), xlim([0 5*ceil(show_horizon/5)])
	plot(0:show_horizon,Cs1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% plot(0:show_horizon,Cs2_show,'-.','linewidth',lw,'color',PLOT.blue);
	plot(0:show_horizon,Cs3_show,'linewidth',lw,'color',PLOT.red);
	plot(0:show_horizon,Cs4_show,':','color',PLOT.purple,'linewidth',lw);
	plot(0:show_horizon,Cs5_show,':','color',PLOT.blue,'linewidth',lw);
	% plot(0:show_horizon,Cs5_show,':','color',PLOT.dgreen,'linewidth',lw);
	title('Savers consumption')
subplot(2,2,2), hold on, grid on, xlabel('Quarters'), ylabel('%'), xlim([0 5*ceil(show_horizon/5)]), ylim([-15 15])
	p1 = plot(0:show_horizon,Y1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% p2 = plot(0:show_horizon,Y2_show,'-.','linewidth',lw,'color',PLOT.blue);
	p3 = plot(0:show_horizon,Y3_show,'linewidth',lw,'color',PLOT.red);
	p4 = plot(0:show_horizon,Y4_show,':','color',PLOT.purple,'linewidth',lw);
	plot(0:show_horizon,Y5_show,':','color',PLOT.blue,'linewidth',lw);
	title('Output')
	% legend([p1 p2 p3], 'no sovereign risk', 'no deleveraging', 'deleveraging', 'location', 'southeast')
	legend([p1 p3], 'no sovereign risk', 'deleveraging', 'location', 'southeast')
subplot(2,2,3), hold on, grid on, xlabel('Quarters'), xlim([0 5*ceil(show_horizon/5)]), ylim([0.85 1.05])
	plot(0:show_horizon,q1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% plot(0:show_horizon,q2_show,'-.','linewidth',lw,'color',PLOT.blue);
	plot(0:show_horizon,q3_show,'linewidth',lw,'color',PLOT.red);
	plot(0:show_horizon,q0_show,'--', 'linewidth',lw,'color',PLOT.dgreen);
	title('Price of government debt')
	legend([p1 p3], 'no sovereign risk', 'deleveraging', 'location', 'southeast')
	% set(gca, 'Ylim', [-14 2])
subplot(2,2,4), hold on, grid on, xlabel('Quarters'), ylabel('% of steady-state GDP'), xlim([0 5*ceil(show_horizon/5)]), ylim([20 34])
	plot(0:show_horizon,T1_show,'--','color',PLOT.yellow,'linewidth',lw);
	% plot(0:show_horizon,T2_show,'-.','linewidth',lw,'color',PLOT.blue);
	plot(0:show_horizon,T3_show,'linewidth',lw,'color',PLOT.red);
	plot(0:show_horizon,T4_show,':','color',PLOT.purple,'linewidth',lw);
	title('Tax Collections')
