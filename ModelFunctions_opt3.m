classdef ModelFunctions_opt3
    methods (Static)
        function R = growth_rate(T,a,T_l,T_m,m)
            %Input =  fitting parameters and temperature T
            %Output =  rate for the given temperature
            
            R = a*T*(T-T_l)*((T_m-T)^(1/m));
        
        end

        function B = birth_rate_suzuki(T,alpha,gamma,lambda,delta,tau)
            %Input =  fitting parameters and temperature T
            %Output =  rate for the given temperature
            
            % Data and function from Ryan et al. (2016)
            Tmin = 5; %Minimum temperature for oviposition
            Tmax = 30; %Maximum temperature for oviposition
            
            if T <= Tmin && T >= Tmax
                B = 0;
            else
                B = alpha * ((gamma + 1)/(pi * lambda^(2 * gamma + 2)) * (lambda^2 - ((T - tau)^2 + delta^2))^gamma);
            end
        end

       function M = mortality_rate(T,a1,b1,c1,d1,e1)
            %Input =  fitting parameters and temperature T
            %Output =  rate for the given temperature
            
            % Data from Ryan et al. (2016) fitted with a fourth polynomial function
            % (script fitmort.py)
        
            if T <=6 && T >= 32 %The function is bounded
                M = 1;
            else
                M= a1 * T^4 + b1 * T^3 + c1 * T^2 + d1 * T + e1;
                
            end
            
            % if M<0
            %    M=0; 
            % end
        end
 

       function Pest_stages = Initialize_stages_ode(Birth,Death,Growth,S_R,R_mate,R_remate,Pest_stages,N_trees,N_stages)
        %Input =  Different rates at each instant
        %Output =  Parameters associated to each stage
        for j=1:N_trees
            for i=1:N_stages %At this point all stages have the same parameters
                Pest_stages(i,j).birth = Birth(j);
                Pest_stages(i,j).death = Death(j);
                Pest_stages(i,j).growth = Growth(j);
                Pest_stages(i,j).sex_ratio = S_R(j);
                Pest_stages(i,j).mate = R_mate(j);
                Pest_stages(i,j).remate = R_remate(j);
              
    
            end
        end

        end



    function [w_d, w_m, w_f] = rate_noise(T,v_a1,v_b1,v_c1,v_d1,v_e1,v_a,v_T_l,v_T_m,v_m,a,T_l,T_m,m)
    
        %Mortality rate error
        w_m= ((v_a1 * T^4)^2 + (v_b1 * T^3)^2 + (v_c1 * T^2)^2 + (v_d1 * T)^2 + (v_e1)^2)^(1/2);
        
        %Growth rate error
        R = a*T*(T-T_l)*((T_m-T)^(1/m));
        
        w_d = (((R/a) * v_a)^2  + ((R/(T-T_l)) * v_T_l)^2 + ( (R/m)*(T_m - T)^(1/(m*(m-1)))  *v_T_m  )^2 + ( (R/m^2) * log(T_m -T)*v_m )^2)^(1/2);
        
        %Fertility rate error (randomly generated based on the other two errors)
        w_f = min(w_m,w_d) + (max(w_m,w_d)-min(w_m,w_d)) .* rand(1,1);
    end

        function trap_list = smart_placement(P_global,N_stages, N_trees,numTrap)
                
                    trace_of_P = zeros(1, N_trees);
                
                    % Compute trace of each per-tree covariance block
                    for i = 1:N_trees
                        idx_start = (i-1)*N_stages + 1;
                        idx_end   = i*N_stages;
                
                        P_block = P_global(idx_start:idx_end, idx_start:idx_end);
                        trace_of_P(i) = trace(P_block);
                    end
                
                    % Sort parcels by descending uncertainty
                    [~, sorted_idx] = sort(trace_of_P, 'descend');
                    trap_list = sorted_idx(1 : numTrap);


               
        end


        
        function [Z,Z_wind,Link]=Adjency_matrix(Columns,Rows)
        %We use as input the number points per column and row
        N=Columns*Rows;
        N = round(N);
        Pos=zeros(N,2);
        pos=1;
        
        %To know file and column
        for i=1:Rows
            for j=1:Columns
                Pos(pos,1) = j;
                Pos(pos,2) =i;
                pos=pos+1; 
            end
        end
        Z=zeros(N,N); %We don't add any extra point
        Z_wind = zeros(N,N);
        Link=zeros(N,2);
        e=1;  
        
        %To stablish which squares ae neighbors and which not
        
        for k=1:N
            for i=1:N
             
                    %To  check if they are neighbours
                    if (Pos(k,1)==Pos(i,1)) 
                        if (abs(Pos(k,2)-Pos(i,2))<2) && (Pos(k,2)~=Pos(i,2))
                            Z(k,i)=1;
                            Link(e,1)=k;
                            Link(e,2)=i;
                            e=e+1;
                            % x coordinate equal
                            if (Pos(k,2)-Pos(i,2)) > 0
                                Z_wind(k,i) = - 1;
                            else
                                Z_wind(k,i) = 1;
                            end
                            
                        end
                    end
                     if (Pos(k,2)==Pos(i,2)) 
                         if (abs(Pos(k,1)-Pos(i,1))<2) && (Pos(k,1)~=Pos(i,1))
                           Z(k,i)=1; 
                           Link(e,1)=k;
                           Link(e,2)=i;
                           e=e+1;
                           
                           % y coordinate equal
                            if (Pos(k,1)-Pos(i,1)) > 0
                                Z_wind(k,i) = - 2;
                            else
                                Z_wind(k,i) = 2;
                            end
                           
                         end
                     end
                
            end
        end


        end

        function [Temp_areas,w] = wind_param(min_t,max_t,Temp_avg,N_trees,wind_speed,wind_direction,saturation,k_factor)
            t_var = min_t + (max_t-min_t) .* rand(length(Temp_avg),N_trees);
            Temp_areas = repmat(Temp_avg,N_trees,1);
            Temp_areas = Temp_areas + t_var';
            check = find(Temp_areas >30);
            Temp_areas(check) = 29.8;
            % Wind direction is degree angles with respect to the north
            wind_x = wind_speed .* cosd(wind_direction);
            wind_y = wind_speed .* sind(wind_direction);
            w=[wind_x', wind_y'];
            % Computation of the wind effect (correcting factor)
            w(find(abs(w)>saturation)) = 0; %If wind more than the saturation limit 0.
            w(isnan(w))=0; %If there is no measurement from the weather station. Wind 0
            w = w*k_factor; %Wind proportional action
        end

    function [A_tot, A_sing] = compute_A_wind_only(N_trees, Pest_stages, stages, Adj_w, w)
    %
    %   POPULATION DYNAMICS + WIND TRANSPORT (NO TRAPS, NO REPELLENT)
  
    
    A_local = [ ...
        -(Pest_stages(1).growth + Pest_stages(1).death),          0,                                  0,                                  0,                                  0,                                  0,                                  0,                            Pest_stages(8).birth;  ... % Egg
    
         Pest_stages(1).growth,                                  -(Pest_stages(2).growth + Pest_stages(2).death), 0,                                  0,                                  0,                                  0,                                  0,                                  0;      ... % L1
    
         0,                                                       Pest_stages(2).growth,             -(Pest_stages(3).growth + Pest_stages(3).death), 0,                                  0,                                  0,                                  0,                                  0;      ... % L2
    
         0,                                                       0,                                  Pest_stages(3).growth,             -(Pest_stages(4).growth + Pest_stages(4).death), 0,                                  0,                                  0,                                  0;      ... % L3
    
         0,                                                       0,                                  0,                                  Pest_stages(4).growth,             -(Pest_stages(5).growth + Pest_stages(5).death), 0,                                  0,                                  0;      ... % Pupa
    
         0,                                                       0,                                  0,                                  0,  (1 - Pest_stages(5).sex_ratio) * Pest_stages(6).growth,   -(Pest_stages(6).growth + Pest_stages(6).death), 0,      0;   ... % Adult male (6)
    
         0,                                                       0,                                  0,                                  0,   Pest_stages(5).sex_ratio * Pest_stages(7).growth,      0, -(Pest_stages(7).mate + Pest_stages(7).death),  0;     ... % NMF (7)  *** FIXED ***
    
         0,                                                       0,                                  0,                                  0,   0,                                                0,   Pest_stages(7).mate, -(Pest_stages(8).growth + Pest_stages(8).death)   ... % MF (8)  *** FIXED ***
    ];
    
    % Replicate into a cell array
    A_sing = repmat({A_local}, 1, N_trees);
    
    %  WIND MOVEMENT MATRIX 
    
    M = zeros(N_trees, N_trees);
    
    if ~isempty(w)
    
        % Horizontal wind (±2)
        if w(1) ~= 0
            idx = (Adj_w == 2  & w(1) > 0) | (Adj_w == -2 & w(1) < 0);
            M(idx) = abs(w(1));
        end
    
        % Vertical wind (±1)
        if numel(w) >= 2 && w(2) ~= 0
            idx = (Adj_w == 1  & w(2) > 0) | (Adj_w == -1 & w(2) < 0);
            M(idx) = abs(w(2));
        end
    
    end
    
    outflow = sum(M,2);   % amount leaving each node
    
    for i = 1:N_trees
        A_sing{i}(6:8, 6:8) = A_sing{i}(6:8, 6:8) - outflow(i) * eye(3);
    end
    
    A_block = blkdiag(A_sing{:});
    

    %  INCOMING WIND

    
    A_jump = zeros(stages);
    A_jump(6,6) = 1;  % AM
    A_jump(7,7) = 1;  % NMF
    A_jump(8,8) = 1;  % MF
    
    A_jump_tot = kron(M.', A_jump);
    
    

    % FINAL MATRIX
    A_tot = A_block + A_jump_tot;
    
    end
    function [x_pred, P_pred] = EKF_predict(x, P, A, Q)
        % Linearized prediction step for EKF
        
        x_pred = A * x;
        P_pred = A * P * A' + Q;
        
        % enforce symmetry
        P_pred = (P_pred + P_pred') / 2;
        end

    function [x_filt, P_filt] = EKF_update(x_pred, P_pred, C, y, R)
    
        S = C * P_pred * C' + R;
        K = P_pred * C' / S;
    
        x_filt = x_pred + K * (y - C * x_pred);
    
        I = eye(size(P_pred));
        P_filt = (I - K*C) * P_pred * (I - K*C)' + K*R*K';   % STABLE
    
        P_filt = (P_filt + P_filt')/2;

    end

    function [whiteIdx, blackIdx, board] = chessboard_grid(rows, cols)
    
        board = zeros(rows, cols);
        whiteIdx = [];
        blackIdx = [];
    
        for i = 1:rows
            for j = 1:cols
                if mod(i + j, 2) == 0
                    board(i, j) = 0;   % white
                    whiteIdx(end+1) = sub2ind([rows, cols], i, j);
                else
                    board(i, j) = 1;   % black
                    blackIdx(end+1) = sub2ind([rows, cols], i, j);
                end
            end
        end
    
%         % Plot
%         figure;
%         imagesc(board);
%         colormap(gray);
%         axis equal tight;
%         grid on;
%         title('Chessboard Grid');
    end



        
    end
end
