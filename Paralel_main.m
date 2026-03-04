clear; clc; close all;

%% PARAMETERS
params = ParametersInputs_opt();
fields = fieldnames(params);
for i = 1:numel(fields)
    eval([fields{i} ' = params.' fields{i} ';']);
end

numTrap = 7;
N_total = N_trees * N_stages;

%% LOAD DATA
data  = readmatrix('Data_rand.xlsx');
data2 = readmatrix('final_x_global.xlsx');
data3 = readmatrix('Data_nonrand.xlsx');
data3 = data3(:,2:end);

Temp_areas = data(:, 5:4:(size(data,2)-2))';
w = data(:, end-1:end);

for k = 2:size(w,1)
    if isnan(w(k,1)), w(k,:) = w(k-1,:); end
end

Simulation_time = size(Temp_areas,2);

%% PRECOMPUTE RATES
Pest_stages(1:N_stages,N_trees) = stage;

G_base = zeros(N_trees, Simulation_time);
B_base = zeros(N_trees, Simulation_time);
M_base = zeros(N_trees, Simulation_time);

for i = 1:N_trees
    for t = 1:Simulation_time
        T = Temp_areas(i,t);
        G_base(i,t) = ModelFunctions_opt3.growth_rate(T,a,T_l,T_m,m);
        B_base(i,t) = ModelFunctions_opt3.birth_rate_suzuki(T,alpha,gamma,lambda,delta,tau);
        M_base(i,t) = ModelFunctions_opt3.mortality_rate(T,a1,b1,c1,d1,e1);
    end
end

%% SPATIAL
[Adj, Adj_w, ~] = ModelFunctions_opt3.Adjency_matrix(sqrt(N_trees), sqrt(N_trees));

%% MONTE CARLO SETTINGS
iter_monte = 100;
x_open_hist_monte = zeros(N_total, Simulation_time, iter_monte);
x_filt_hist_monte = zeros(N_total, Simulation_time, iter_monte);
traceP_monte      = zeros(Simulation_time, iter_monte);
traceP_adult_monte= zeros(Simulation_time, iter_monte);

Q_global = eye(N_total)*0.1;

 parpool(); 
tic,
parfor monte = 1:iter_monte


    x_global = zeros(N_total,1);
    P_global = eye(N_total);
    Q_global = 0.1 * eye(N_total);

    traceP_local       = zeros(Simulation_time,1);
    traceP_adult_local = zeros(Simulation_time,1);

    for i = 1:N_trees
        x_i = zeros(N_stages,1);
        stage_idx = randi(N_stages);
        x_i(stage_idx) = randi([80 100]);
        x_global((i-1)*N_stages+1:i*N_stages) = x_i;
    end

    x_open_hist = zeros(N_total, Simulation_time);
    x_filt_hist = zeros(N_total, Simulation_time);

    x_open_hist(:,1) = x_global;
    x_filt_hist(:,1) = x_global;

    x_pred = x_global;

    % INITIAL A 
    Pest_stages_local = ModelFunctions_opt3.Initialize_stages_ode( ...
        B_base(:,1), M_base(:,1), G_base(:,1), ...
        s_r*ones(N_trees,1), r_mate*ones(N_trees,1), r_remate*ones(N_trees,1), ...
        Pest_stages, N_trees, N_stages);

    A_cont = ModelFunctions_opt3.compute_A_wind_only( ...
        N_trees, Pest_stages_local, N_stages, Adj_w, w(1,:));

    A_dis = c2d(ss(A_cont,[],eye(N_total),[]),1).A;
      %% Grid chess
%       [whiteIdx, blackIdx, board] = ModelFunctions_opt3.chessboard_grid(sqrt(N_trees), sqrt(N_trees));


    % Simulation loop starts
    for t = 1:Simulation_time-1

        % Open loop
        x_global = A_dis * x_global;
        x_open_hist(:,t+1) = x_global;

        % EKF predict
        [x_pred, P_pred] = ModelFunctions_opt3.EKF_predict( ...
            x_pred, P_global, A_dis, Q_global);

        x_filt_hist(:,t+1) = x_pred;
        P_global = P_pred;
         

        % Measurement update
        if mod(t,7) == 0
         %% Smart Startegy
%         observed = ModelFunctions_opt3.smart_placement(P_global,N_stages, N_trees,numTrap);
        %% Random Startegy
        observed = randperm(N_trees, numTrap);
        %% Chess Straegy
%           sensor = [4 4; 2 2; 6 2; 2 6; 6 6; 1 4; 7 4];
%           observed = sub2ind([sqrt(N_trees) sqrt(N_trees)], sensor(:,1), sensor(:,2));
        %% Diagonal and Anti-Diagonal strategy
%           Diagonal = [25     9    13    37    41    22    28]; % Diag
%           AntiDiagonal = [1     9    17    25    33    41    49]; %Anidiag
%           observed  = Diagonal;

            m_meas = numTrap * 3;
            C = zeros(m_meas, N_total);
            y = zeros(m_meas,1);

            row = 1;
            for k = 1:numTrap
                tree = observed(k);
                col = (tree-1)*N_stages + 6;
                C(row:row+2, col:col+2) = eye(3);

                dcol = (tree-1)*4 + 1;
                y(row:row+2) = data3(t+1, dcol:dcol+2)';
                row = row + 3;
            end
            % CORRECTION
            [x_filt, P_filt] = ModelFunctions_opt3.EKF_update( ...
                x_filt_hist(:,t+1), P_global, C, y, 0.1*eye(m_meas));

            x_filt_hist(:,t+1) = x_filt;
            x_pred = x_filt;
            P_global = P_filt;
        end

        % TRACE 
        adult_idx = reshape((1:N_trees)'*N_stages + (6:8) - N_stages,[],1);
        P_adult = P_global(adult_idx, adult_idx);

        traceP_adult_local(t+1) = trace(P_adult);
        traceP_local(t+1) = trace(P_global);

        %  UPDATE A 
        Pest_stages_local = ModelFunctions_opt3.Initialize_stages_ode( ...
            B_base(:,t+1), M_base(:,t+1), G_base(:,t+1), ...
            s_r*ones(N_trees,1), r_mate*ones(N_trees,1), r_remate*ones(N_trees,1), ...
            Pest_stages_local, N_trees, N_stages);

        A_cont = ModelFunctions_opt3.compute_A_wind_only( ...
            N_trees, Pest_stages_local, N_stages, Adj_w, w(t+1,:));

        A_dis = c2d(ss(A_cont,[],eye(N_total),[]),1).A;
    end

    % STORE 
    x_open_hist_monte(:,:,monte) = x_open_hist;
    x_filt_hist_monte(:,:,monte) = x_filt_hist;
    traceP_monte(:,monte)        = traceP_local;
    traceP_adult_monte(:,monte)  = traceP_adult_local;
end
elapsed_time = toc;
disp(elapsed_time)  

%% SAVE 
save('MonteCarloRand.mat', ...
     'x_open_hist_monte', 'x_filt_hist_monte', ...
     'traceP_monte', 'traceP_adult_monte');

%% PLOTS
figure; hold on; grid on;

plot(mean(traceP_monte,2),        'k', 'LineWidth', 2);
plot(mean(traceP_adult_monte,2),  'r', 'LineWidth', 2);
figure; hold on; grid on;

plot(mean(traceP_monte,2),        'k', 'LineWidth', 2);
plot(mean(traceP_adult_monte,2),  'r', 'LineWidth', 2);


all_values = traceP_adult_monte(:);   

figure; grid on
boxplot(all_values, 'Notch','on');
grid on;

ylabel('Trace(P_{adult})');
title('Monte Carlo Distribution of EKF Covariance (All times)');


% Mean over Monte Carlo
adult_est_mean  = mean(x_filt_hist_monte, 3);
adult_open_mean = mean(x_open_hist_monte, 3);

%%
meas_idx = find(mod(1:Simulation_time,7) == 0);

for parcel = 1:49
    figure; hold on; grid on;

    % EKF + Open-loop (sum of stages 6–8) 
    idx_state = (parcel-1)*N_stages + (6:8);

    adult_est_sum  = sum(adult_est_mean(idx_state,:), 1);
    adult_open_sum = sum(adult_open_mean(idx_state,:), 1);

    plot(adult_est_sum,  'r', 'LineWidth', 2);
    plot(adult_open_sum, 'b--','LineWidth', 2);

    %  Measurements (data3: 8 cols per parcel, adults = 6:8) 
    cols_adults = (parcel-1)*4+ (1:3);
    adult_meas_sum = sum(data3(meas_idx, cols_adults), 2);

    scatter(meas_idx, adult_meas_sum, 40, 'go', 'filled');

    legend('EKF mean','Open-loop mean','Measurements');
    xlabel('Time');
    ylabel('Adults (sum stages 6–8)');
    title(['Parcel ', num2str(parcel)]);
end
%%
load('MonteCarloRand.mat',      'traceP_monte');
trace_rand = traceP_monte;

load('MonteCarloSmart.mat',     'traceP_monte');
trace_smart = traceP_monte;

load('MonteCarloChess.mat',     'traceP_monte');
trace_chess = traceP_monte;

load('MonteCarloDiagonal.mat',  'traceP_monte');
trace_diag = traceP_monte;

% final_rand  = log10(trace_rand(end,:));
% final_smart = log10(trace_smart(end,:));
% final_chess = log10(trace_chess(end,:));
% final_diag  = log10(trace_diag(end,:));

data = [a(:), ...
        b(:), ...
        c(:), ...
        d(:)];

figure; boxplot(data, ...
    'Labels', {'Random','Smart','Chessboard','Diagonal'}, ...
    'Notch','on');

grid on;
ylabel('log_{10}( Trace(P_{adult}) )');
title('Monte Carlo EKF Uncertainty Comparison (Final Time)');

%%
delete(gcp);
