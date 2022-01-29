%% load 2D dataset
%load('Waves_dp_100.mat');
load('Waves_dp.mat');
%% load 1D dataset
%load('Burgers_idp.mat');
%load('Heat_dp.mat');
%% Build normal dictionaries
%OneD_dict;
TwoD_dict;
%% Build transformed dicitonaries, for 2D waves equation dataset in Sec. III-D
% add_noise;
% TwoD_dict_int;  % use intergration transformation  
%% Recovery
lasso_seq; % if using transformed dictionary, make sure the input to lasso_seq is changed to "Integral_res"
recover_c_a_2dwave; % not used for 1D cases


