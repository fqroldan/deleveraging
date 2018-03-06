%% This file produces pictures detailing parts of the savers welfare
% Francisco Roldán, NYU. October 2016
clear all
close all
clc
%%% Remember to set run_outside = 1 in comparepolicies.m before running %%%

%% First run comparepolicies without private deleveraging
withpriv = 0;
save('withpriv.mat', 'withpriv')
comparepolicies

% Save results
Cs_without_priv = Cs_mat(2,1:show_horizon)';
Cs_without_priv_def = Csd_mat(2,1:show_horizon)';

Ns_without_priv = Ns_mat(2,1:show_horizon)';
Ns_without_priv_def = Nsd_mat(2,1:show_horizon)';

% Store results
save('temp.mat')

clear all

%% Now run comparepolicies again, with private deleveraging
withpriv = 1;
save('withpriv.mat', 'withpriv')
comparepolicies

% Save results
Cs_with_priv = Cs_mat(2,1:show_horizon)';
Cs_with_priv_def = Csd_mat(2,1:show_horizon)';

Ns_with_priv = Ns_mat(2,1:show_horizon)';
Ns_with_priv_def = Nsd_mat(2,1:show_horizon)';

% Load results from before
load temp.mat
close all
%% Plot together
figure('position',[476   385   886   570])
subplot(121), hold on, grid on, xlabel('Quarters')
plot(Cs_without_priv,'linewidth',1.5), plot(Cs_with_priv,'linewidth',1.5)
legend('With priv del', 'Without priv del', 'location','northeast')
subplot(121), hold on, plot(Cs_without_priv_def,':','color',PLOT.blue,'linewidth',1.5),plot(Cs_with_priv_def,':','color',PLOT.red,'linewidth',1.5)
title('Consumption - Savers')

subplot(122), hold on, grid on, xlabel('Quarters')
plot(Ns_without_priv,'linewidth',1.5), plot(Ns_with_priv,'linewidth',1.5)
legend('Without priv del', 'With priv del', 'location','southeast')
subplot(122), hold on, plot(Ns_without_priv_def,':','color',PLOT.blue,'linewidth',1.5),plot(Ns_with_priv_def,':','color',PLOT.red,'linewidth',1.5)
title('Labor supply - Savers')

suptitle('Govt deleveraging case')

% eval(['print -dpng ../../Graphs/Nperiods/savers_comparison.png -r200'])
