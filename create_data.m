% x_open_hist : 392 × T  (49 trees × 8 stages)
% w           : T × 2
% We build a T × (49×4) matrix
% This code was developped by Benhamouche Ouassim and Luca Rossini
%% DATA3
T = size(x_open_hist, 2);
numTrees = 49;
numStages = 8;
% Temp_areas = Temp_areas';

FinalAdultsWind = zeros(T, numTrees*4);
for i = 1:numTrees
    temp_col = 2 + (i - 1) * 4 + 3; 
    Temp_areas(:,i) = a(:,temp_col);
end

for tree = 1:numTrees
    
    % Adult indices in the global vector
    idx_AM  = (tree-1)*numStages + 6;
    idx_NMF = (tree-1)*numStages + 7;
    idx_MF  = (tree-1)*numStages + 8;
    
    % Column group index in output matrix
    out_col = (tree-1)*4 + 1;
    
    % Fill output
    FinalAdultsWind(:, out_col    ) = x_open_hist(idx_AM , :)';   % AM
    FinalAdultsWind(:, out_col + 1) = x_open_hist(idx_NMF, :)';   % NMF
    FinalAdultsWind(:, out_col + 2) = x_open_hist(idx_MF , :)';   % MF
    FinalAdultsWind(:, out_col + 3) = Temp_areas(:,tree);                      % wind column 1
end
FinalAdultsWind = real(FinalAdultsWind);


a(186:197,198:199) = a(174:185,198:199);
a=real(a);
FinalAdultsWind(:,end+1:end+2) = a(:,198:199);
save('data','FinalAdultsWind')