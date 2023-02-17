% This script provides a demonstration of the Bouc-Wen class model proposed by Oh et al. (2023).
% Experimental dataset is adopted from the Structural Performance Database complied by the Pacific Earthquake Engineering Research (PEER) center
% https://nisee.berkeley.edu/spd/ (Accessibility confirmed by 02/16/23)
% First dataset is used; Gill, W. D., Park, R., & Priestley, M.J.N. Ductility of Rectangular Reinforced Concrete Columns With Axial Load. Report 79-1, Department of Civil Engineering, University of Canterbury, Christchurch, New Zealand, February 1979, 136 pages.

%% clearing
clear; clc; close all;

%% load the experimental data
data = load('experimental data.mat');
data = data.data_exp;

disp  = data(:,1);
force = data(:,2);

%% load parameter for the hmBWBN model
params = load('params.mat');
params = params.params_est;

%% compare the experimental data and the model

force_est = BoucWen(params,disp);

figure;
plot(disp,force,'k-','linewidth',1.2); grid on; hold on;
plot(disp,force_est,'r--','linewidth',1.2);
xlabel('Displacement (m)'); ylabel('Lateral load (g)');
legend('Experimental data','hmBWBN model','location','Southeast');
set(gca,'fontname','Times New Roman','fontsize',13);
