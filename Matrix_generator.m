clear

% loading x_global_hist & matrix with temp and wind 

histy2 = load('x_global_hist_m4.mat');

% histy2 = table2array(histy2);
histy2 = histy2.x_open_hist';
histy2 = real(histy2);
data = load('data.mat');
data = data.FinalAdultsWind;
data = real(data);
% Data_Matrix_final_rand = load('newdatawass4.mat');
% Data_Matrix_final_rand = Data_Matrix_final_rand.c;%table2array(Data_Matrix_final_rand);
% 
% Selecting the adult stages + temperature column;
num_nodes = 49; 
num_stages = 8; 

final_matrix3 = nan(size(histy2, 1), num_nodes * 4); % 4 columns per node: 6, 7, 8 + empty.

for node_idx = 1:num_nodes

    base_idx = (node_idx - 1) * num_stages; 
    col_6 = base_idx + 6; 
    col_7 = base_idx + 7; 
    col_8 = base_idx + 8; 

    final_col_start = (node_idx - 1) * 4 + 1; 
    final_matrix3(:, final_col_start) = histy2(:, col_6); 
    final_matrix3(:, final_col_start + 1) = histy2(:, col_7); 
    final_matrix3(:, final_col_start + 2) = histy2(:, col_8); 
    final_col_start + 3; 
end

% Randomizing the matrix (uncomment only when we want to generate a randomized matrix)

random_factors = 0.2;

final_matrix3 = final_matrix3 + final_matrix3 * unifrnd(-0.2,0.2);

% Final matrix (with time-steps, non-rand temp and wind)(data could be one of those matrix i did, like Data_rand)

[n_rows, n_cols] = size(final_matrix3); 
n_nodes = round(n_cols / 4);               

%Time-steps in the first column
time_steps = (1:n_rows)';             
Final_Matrix = zeros(n_rows, n_cols + 3); 
Final_Matrix(:, 1) = time_steps; 

% Adult stages addition
Final_Matrix(:, 2:n_cols+1) = final_matrix3;

%Temperatures addition 
for i = 1:n_nodes
    temp_col = 2 + (i - 1) * 4 + 3; 
    Final_Matrix(:, temp_col) = data(:,temp_col-1);
end

% Wind direction and wind speed addition 
Final_Matrix(:, end-1:end) = data(:, end-1:end); 
% Final_Matrix = Final_Matrix(:,1:end-1);


% writematrix(x_global_hist, 'x_global_hist.xlsx');
writematrix(histy2, 'final_x_global.xlsx');
writematrix(Final_Matrix, 'Data_nonrand.xlsx')
writematrix(Final_Matrix, 'Data_rand.xlsx')  %in un ciclo attivo quella per stampare la rand e nell'altro stampo la non rand)

