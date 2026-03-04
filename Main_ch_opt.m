% ODE model for Drosophila suzukii
% Time varying (based on temperature) – OPEN LOOP ONLY
clear; clc; close all;
% rng(1)

%  PARAMETERS

params = ParametersInputs_opt();
fields = fieldnames(params);
for i = 1:numel(fields)
    eval([fields{i} ' = params.' fields{i} ';']);
end
numTrap = 7;
N_total = N_trees * N_stages;

%  LOAD DATA
noise = 0.2;
data  = readmatrix('Data_rand.xlsx');
data2 = readmatrix('final_x_global.xlsx');   % ground truth
data3 = readmatrix('Data_nonrand.xlsx');     % NOT used here, only EKF

data3 = data3(:,2:end);

Temp_areas = data(:, 5:4:(size(data,2)-2))'; % temperatures per tree x time
w = data(:, end-1:end);                      % wind



for row = 2:size(w,1)
    if isnan(w(row,1)), w(row,:) = w(row-1,:); end
end

% Extract adults ground truth (for plotting later)
startCols = 2:4:(size(data2, 2) - 3);
groupIndices = [];
for col = startCols
    groupIndices = [groupIndices, col:(col+2)];
end

Simulation_time = size(Temp_areas,2);

%% PRECOMPUTE RATES (growth, birth, death) FOR EACH TREE/TIME
%  
Pest_stages(1:N_stages,N_trees) = stage; % We create an array of the stage class

G_base = zeros(N_trees, Simulation_time);
B_base = zeros(N_trees, Simulation_time);
M_base = zeros(N_trees, Simulation_time);

for i = 1:N_trees
    for t = 1:Simulation_time
        T = Temp_areas(i,t);
        G_base(i,t) = ModelFunctions_opt3.growth_rate( ...
            T, a, T_l, T_m, m);
        B_base(i,t) = ModelFunctions_opt3.birth_rate_suzuki( ...
            T, alpha, gamma, lambda, delta, tau);
        M_base(i,t) = ModelFunctions_opt3.mortality_rate( ...
            T, a1, b1, c1, d1, e1);
    end
end



%% SPATIAL CONFIGURATION

[Adj, Adj_w, Link] = ModelFunctions_opt3.Adjency_matrix( ...
    sqrt(N_trees), sqrt(N_trees));

iter_monte = 1;
x_global_hist_monte = zeros(N_total, Simulation_time, iter_monte);
P_global = eye(N_total)*1;
Q_global = eye(N_total)*0.1;
tic;

for monte = 1:iter_monte
    x_global = [];
    for i = 1:N_trees
        x_i = zeros(N_stages,1);
        stage_idx = randi([1 N_stages]);
        x_i(stage_idx) = randi([80 100]);   % initial individuals
        x_global = [x_global; x_i];
    end

    x_open_hist = zeros(N_total, Simulation_time);
    x_open_hist(:,1) = data2(1,:)';
    x_pred = data2(1,:)';
    x_filt_hist = x_open_hist;
%     x_global = data2(1,:)';

    %% ---------- INITIAL A (continuous & discrete) ----------
    Pest_stages = ModelFunctions_opt3.Initialize_stages_ode( ...
        B_base(:,1), M_base(:,1), G_base(:,1), ...
        s_r*ones(N_trees,1), ...
        r_mate*ones(N_trees,1), ...
        r_remate*ones(N_trees,1), ...
        Pest_stages, N_trees, N_stages);


   A_cont = ModelFunctions_opt3.compute_A_wind_only( ...
                 N_trees, Pest_stages, N_stages, Adj_w, w(1,:));


    %  STABILITY CHECK AT t = 1
    lam = eig(A_cont);
    fprintf('max Re(eig(A_cont)) at t=1: %.4f\n', max(real(lam)));

    sysc = ss(A_cont, [], eye(N_total), []);
    sysd = c2d(sysc, 1, 'zoh');    % dt = 1
    A_dis_global = sysd.A;
    %% Grid chess
    [whiteIdx, blackIdx, board] = ModelFunctions_opt3.chessboard_grid(sqrt(N_trees), sqrt(N_trees));

    %% Simulation
    for t = 1:Simulation_time-1
        x_global = A_dis_global * x_global;
        x_open_hist(:, t+1) = x_global;        
    %% EKF    
    % 1) PREDICTION   
    [x_pred, P_pred] = ModelFunctions_opt3.EKF_predict(x_pred, P_global, A_dis_global, Q_global);
    x_filt_hist(:,t+1) = x_pred;                   % EKF prediction storage
    P_global = P_pred;


    
    if mod(t,7) == false
        % 2) MEASUREMENT MODEL 
        numTrees_obs = numTrap;
        numObsPerTree = 3;       % AM, NMF, MF
        m_meas = numTrees_obs * numObsPerTree;
        %% Smart Startegy
%         observed = ModelFunctions_opt3.smart_placement(P_global,N_stages, N_trees,numTrap);
        %% Random Startegy
        observed = randperm(N_trees, numTrap);
        %% Chess Straegy
%           sensor = [4 4; 2 2; 6 2; 2 6; 6 6; 1 4; 7 4];
%           observed = sub2ind([sqrt(N_trees) sqrt(N_trees)], sensor(:,1), sensor(:,2));
        %% Diagonal and Anti-Diagonal strategy
%           Diagonal = [25     9    13    37    41    22    28]; % Diag
%           observed  = Diagonal;
%           AntiDiagonal = [1     9    17    25    33    41    49]; %Anidiag
%           observed  = AntiDiagonal;
         %% Observations
        C_global = zeros(m_meas, N_total);
        y_vec = zeros(m_meas, 1);
        
        row = 1;
        
        for k = 1:numTrees_obs
            tree = observed(k);
        
            % Stages 6–8 of this tree
            col_start = (tree-1)*N_stages + 6;  % AM, NMF, MF in state vector
            C_global(row:row+2, col_start:col_start+2) = eye(3);

        
            % Get adult data from file
            adult_col_start = (tree-1)*4 + 1;
            y = data3(t+1, adult_col_start : adult_col_start + 2)';
       
            y_vec(row:row+2) = y;
            row = row + 3;

        end

        % 3) UPDATE 
    
        R_global = 0.1 * eye(m_meas);
    
        [x_filt, P_filt] = ModelFunctions_opt3.EKF_update(x_filt_hist(:,t+1), P_global, C_global, y_vec, R_global);
    
        x_filt_hist(:,t+1) = x_filt;
        x_pred = x_filt; 
        P_pred = P_filt;
        % Replace prior with posterior
        P_global = P_filt;


    end

    adult_idx = [];
    for i = 1:N_trees
        adult_idx = [adult_idx, (i-1)*N_stages + (6:8)];
    end
    
    P_adult = P_global(adult_idx, adult_idx);
    traceP_adult(t+1) = sum(diag(P_adult));   % sum of variances of all adult states
    traceP(t+1) = (trace(P_adult));


        % UPDATE PARAMETERS 
        Pest_stages = ModelFunctions_opt3.Initialize_stages_ode( ...
            B_base(:,t+1), M_base(:,t+1), G_base(:,t+1), ...
            s_r*ones(N_trees,1), ...
            r_mate*ones(N_trees,1), ...
            r_remate*ones(N_trees,1), ...
            Pest_stages, N_trees, N_stages);



        A_cont = ModelFunctions_opt3.compute_A_wind_only( ...
                     N_trees, Pest_stages, N_stages, Adj_w, w(t+1,:));

        sysc = ss(A_cont, [], eye(N_total), []);
        sysd = c2d(sysc,1, 'zoh');
        A_dis_global = sysd.A;

    
        
    end
    x_global_hist_monte(:,:,monte) = x_open_hist;
end

toc;

%%
adult_true = data2(:, groupIndices)';   % ground truth adults

adult_open = [];

for tree = 1:N_trees
    idx_start = (tree-1)*N_stages + 6;   % stages 6–8 = adults
    idx_end   = idx_start + 2;
    adult_open = [adult_open; x_open_hist(idx_start:idx_end, :)];
end

% for tree = 1:numTrees
%     idx_start_true = (tree-1)*3 + 1;   % AM, NMF, MF positions in adult_true
% 
%     % In data3 structure: adults begin at column 2+(tree-1)*4
%     idx_data3 = 2 + (tree-1)*4;
% 
%     adult_meas(idx_start_true : idx_start_true+2, :) = data3(:, idx_data3:idx_data3+2)';
% end
% adult_true = adult_meas

adult_est  = [];

for tree = 1:N_trees
    idx_start = (tree-1)*N_stages + 6;
    idx_end   = idx_start + 2;
    adult_est  = [adult_est;  x_filt_hist(idx_start:idx_end, :)];  
end

adult_meas = zeros(3*N_trees, Simulation_time);
for tree = 1:N_trees
    idx_true  = (tree-1)*3 + 1;      % in adult_true
    idx_data3 =  1 + (tree-1)*4;      % in data3
    adult_meas(idx_true:idx_true+2, :) = data3(:, idx_data3:idx_data3+2)';
end


for parcel = 5 : 20   % change this to inspect another parcel (1..49)

    figure; hold on; grid on
    
    corr_days = find(mod(1:Simulation_time,7)==0);
    
    parcel_idx = (parcel-1)*3 + 3;   % MF only
    
    true_vals = adult_true(parcel_idx, corr_days);
    meas_vals = adult_meas(parcel_idx, corr_days);
    
    % NORMAL LINES
    plot(adult_true(parcel_idx,:), 'k', 'LineWidth', 2);
    plot(adult_est(parcel_idx,:),  'r', 'LineWidth', 2);
    plot(adult_open(parcel_idx,:), 'b--', 'LineWidth', 2);
    
    % SCATTERS
    scatter(corr_days, true_vals, 60, 'ko', 'filled');         % true
    scatter(corr_days, meas_vals, 60, 'go', 'filled');         % measured

end

%%

figure; hold on;
plot(traceP, 'LineWidth', 2);
xlabel('Time');
ylabel('Trace(P)');
title('Evolution of EKF Covariance (trace)');
grid on;


