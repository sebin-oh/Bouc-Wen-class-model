% Created: April 14, 2025
% Author : Sebin Oh @ University of California, Berkeley
% Contact: sebin.oh@berkeley.edu
% Description: 
% This script estimates parameters of the Bouc-Wen class models for a given hysteresis curve using genetic algorithm.
% While the Bouc-Wen model proposed in Oh et al.(2023) is adopted here for illustrative purposes, any other model is applicable.

%%
clear; clc; close all;

%% Load hysteresis curve
data = load("experimental data.mat");
data = data.data_exp;

disp = data(:,1);
force = data(:,2);

% figure;
% plot(disp*1e2, force, 'k-', 'linewidth', 1.2); grid on;
% xlabel('Displacement (cm)'); ylabel('Load (g)');
% set(gca,'fontname','Times New Roman','fontsize',11);

%% Genetic Algorithm
para_estimate = 15; % Number of parameters to estimate

% Basic parameters for GA
num_population = 500;
num_parents = ceil(num_population*0.2);
num_GA = 100;
LR = -log(0.01)/num_GA; % make the weight of the mutation goes to 0.1 at terminal value

% Load T, F_y, and alpha values estimated from a pushover curve
initparam = load("initparam.mat");
initparam = initparam.initparam;

% Set bounds for the parameters
eps = 1e-5;
bounds = [initparam(1)*(1+eps),initparam(1)*(1-eps); % 3, 0.05; % 1: period
    initparam(2)*(1+eps),initparam(2)*(1-eps); % 1, 0.05; % 2: F_y
    initparam(3)*(1+eps),initparam(3)*(1-eps); % 0.5, 0; % 3: alpha
    1,0.01; % 4: beta
    5,1; % 5: n
    1.0,0; % 6: del_nu
    1.0,0; % 7: del_eta
    0.9999,0.0; % 8: zeta_0
    10,0; % 9: p
    0.3,0.01; % 10: q
    10.0,0.1; % 11: psi
    3.0,0.0; % 12: del_psi
    10.0,0.01; % 13: lambda
    200,0; % 14: c_eps
    3.0,0.05]; % 15: c_pch

% Initial population
pop = transpose(bounds(:,1)-bounds(:,2)).*rand(num_population, para_estimate) ...
    +transpose(bounds(:,2));

% Iterations
for jj =1:num_GA

    fprintf('%dth generation\n',jj);

    fitness = GA_fitness(pop, disp, force,bounds); % Calculate objective function values

    parents = GA_mating(pop, fitness, num_parents, bounds); % Find the candidates for parents --> the fitness values are smallest

    offspring = GA_crossover(parents, num_population-num_parents); % perform crossover

    crossover = GA_mutation(offspring, LR, jj, bounds); % perform mutation

    pop = [parents; crossover];

end

% Find the best chromosomes
fitness = GA_fitness(pop, disp, force,bounds); % calculate fitness values
[~,Best_idx] = min(fitness);
Best_set = pop(Best_idx,:);

%%
force_est = BoucWen(Best_set,disp);

figure;
plot(disp*1e2, force, 'k-', 'linewidth', 1.2); grid on; hold on;
plot(disp*1e2, force_est,'r--','linewidth',1.2);
xlabel('Displacement (cm)'); ylabel('Load (g)');
legend('Experimental data','Fitted BW model','location','Southeast');
set(gca,'fontname','Times New Roman','fontsize',11);
