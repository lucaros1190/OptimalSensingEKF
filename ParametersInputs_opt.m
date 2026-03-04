classdef ParametersInputs_opt
    properties
      N = 10000;
    sampling_time = 1;
    rep_efficiency = 0;
    trap_efficiency = 0;
    mort_trap = 0;
    % Wind
    k_factor = 1/500; % Proportional factor wind
    saturation = 5; % Wind at which there is no movement

    v_a1 = 1.4E-06;
    v_b1 = 0.0008;
    v_c1 = 0.02;
    v_d1 = 0.3E-05;
    v_e1 = 0.9;

    %Uncertainty growth rate
    v_a = 0.15 * (10^(-4));
    v_T_l = 2;
    v_T_m = 1;
    v_m = 3;
    
     %Growth rate (Drosophila suzukii)
    a = 1.2*(10^(-4));
    T_l = 3;
    T_m = 30;
    m = 6;
    %  
    %Mortality rate (Drosophila suzukii)
    a1 = -2.98E-06;    %-5.4E-06;
    b1 =0.0004020;   %0.0005194;
    c1 = -0.0104770;  %-0.0116827;
    d1 = 1.70E-05;   %2.16E-05;
    e1 = 1.2080000;  %1.3146586;


    % Birth rate (Drosophila suzukii)
    alpha = 659.06;
    gamma = 88.53;
    lambda = 52.32;
    delta = 6.06;
    tau = 22.87;
 % Sex ratio
    s_r = 0.5;
    % Mating ratio
    r_remate = 0;
    r_mate = 1;
    % Additional parameters
    N_stages = 8; % Eggs/3 larva stages/ pupa /male/ unmated female/ mated female + traps(male/ unmated female/ mated female)
    N_trees = 49; %Number of tress
    max_t = 0.2;
    min_t = -0.2;
    counter = 1;
    counter_meas = 1;
    counter_check = 1;
    trap_on = 0;
    w_d=0.2; % We put them initially to 0
    w_m=0.2;
    w_f=0.2;

    trap_time = 8;
    trap_mort= 0.0;
    trap_attrac = 0;
    noobs_length = 1;
    end
end