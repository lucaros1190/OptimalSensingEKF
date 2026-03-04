clear; clc; close all;
rng(1)

%% select the name of thepolicy
policy = "smart";   
% options: "smart", "random", "chess", "diagonal", 
numTrap_list = [5 15 20];          
iter_monte   = 1;   

%% PARAMETERS
params = ParametersInputs_opt();
fields = fieldnames(params);
for i = 1:numel(fields)
    eval([fields{i} ' = params.' fields{i} ';']);
end

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

%% SETTINGS
Q_global = 0.1 * eye(N_total);

traceP_all       = zeros(length(numTrap_list), Simulation_time);
traceP_adult_all = zeros(length(numTrap_list), Simulation_time);

% parpool();

for nt = 1:length(numTrap_list)

    numTrap = numTrap_list(nt);

    x_global = zeros(N_total,1);
    P_global = eye(N_total);

    for i = 1:N_trees
        x_i = zeros(N_stages,1);
        stage_idx = randi(N_stages);
        x_i(stage_idx) = randi([80 100]);
        x_global((i-1)*N_stages+1:i*N_stages) = x_i;
    end

    x_pred = x_global;

    traceP_local       = zeros(Simulation_time,1);
    traceP_adult_local = zeros(Simulation_time,1);

    Pest_stages_local = ModelFunctions_opt3.Initialize_stages_ode( ...
        B_base(:,1), M_base(:,1), G_base(:,1), ...
        s_r*ones(N_trees,1), r_mate*ones(N_trees,1), r_remate*ones(N_trees,1), ...
        Pest_stages, N_trees, N_stages);

    A_cont = ModelFunctions_opt3.compute_A_wind_only( ...
        N_trees, Pest_stages_local, N_stages, Adj_w, w(1,:));

    A_dis = c2d(ss(A_cont,[],eye(N_total),[]),1).A;

    for t = 1:Simulation_time-1

        x_global = A_dis * x_global;

        [x_pred, P_pred] = ModelFunctions_opt3.EKF_predict( ...
            x_pred, P_global, A_dis, Q_global);

        P_global = P_pred;

        if mod(t,7) == 0

        switch policy
        
            case "smart"
                observed = ModelFunctions_opt3.smart_placement( ...
                    P_global, N_stages, N_trees, numTrap);
        
            case "random"
                observed = randperm(N_trees, numTrap);
        
            case "chess"
                sensor = [4 4; 2 2; 6 2; 2 6; 6 6; 1 4; 7 4];
                observed = sub2ind([sqrt(N_trees) sqrt(N_trees)], ...
                                   sensor(:,1), sensor(:,2));
        
            case "diagonal"
                observed = [25 9 13 37 41 22 28];
        
            case "antidiagonal"
                observed = [1 9 17 25 33 41 49];
        
            otherwise
                error("Unknown sensing policy");
        end


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

            [x_filt, P_filt] = ModelFunctions_opt3.EKF_update( ...
                x_pred, P_global, C, y, 0.1*eye(m_meas));

            x_pred  = x_filt;
            P_global = P_filt;
        end

        adult_idx = reshape((1:N_trees)'*N_stages + (6:8) - N_stages,[],1);
        P_adult = P_global(adult_idx, adult_idx);

        traceP_local(t+1)       = trace(P_global);
        traceP_adult_local(t+1) = trace(P_adult);

        Pest_stages_local = ModelFunctions_opt3.Initialize_stages_ode( ...
            B_base(:,t+1), M_base(:,t+1), G_base(:,t+1), ...
            s_r*ones(N_trees,1), r_mate*ones(N_trees,1), r_remate*ones(N_trees,1), ...
            Pest_stages_local, N_trees, N_stages);

        A_cont = ModelFunctions_opt3.compute_A_wind_only( ...
            N_trees, Pest_stages_local, N_stages, Adj_w, w(t+1,:));

        A_dis = c2d(ss(A_cont,[],eye(N_total),[]),1).A;
    end

    traceP_all(nt,:)       = traceP_local;
    traceP_adult_all(nt,:) = traceP_adult_local;
end

delete(gcp);

%% PLOTS — OVERALL TRACE VS NUMBER OF TRAPS
figure; hold on; grid on;
plot(numTrap_list, mean(traceP_all,2), 'k-o','LineWidth',2);
plot(numTrap_list, mean(traceP_adult_all,2), 'r-o','LineWidth',2);
xlabel('Number of traps');
ylabel('Mean Trace');
legend('Trace(P)','Trace(P_{adult})');
title('Overall EKF Uncertainty vs Number of Traps');

%%
filename = sprintf('MonteCarlo_%s.mat', policy);

save(filename, ...
     'x_open_hist_monte', ...
     'x_filt_hist_monte', ...
     'traceP_monte', ...
     'traceP_adult_monte');


