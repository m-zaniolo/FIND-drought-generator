function param = generate_param_structure()
% this function contains a number of parameters related to drought
% definition and optimization. Optimization parameters may require manual
% tuning to ensure optimization convergence and speed. Refer to the paper
% for guidelined on how to set these parameters.

%% parameters used for experiment 1 and 2 
param.min_drought_duration  = 32; % minimum drought duration
param.min_drought_intensity = -0.6; % minimum drought intensity
param.ssi_time_scale        = 12; % aggregation time in months for SSI calculation 
param.nmonth                = 48; % portion of the time series that is initially replaced in months
param.nmonths_end_drought   = 12; % number of consecutive non-negative SSI values to end a drought

param.nyear_togenerate     = 100; % length of the generated scenarios in years
param.H                    = param.nyear_togenerate*12; % % length of the generated scenarios in months

% define annealing schedule. 
param.Nm = 600; % number of swaps for each temperature and nmonth
param.m  = 15 ; % number of temperature and nmonth reductions
param.T  = 0.0001; % initial temperature
param.decrease_rate = 0.8; %decrease rate of temperature and nmonth

%define weights 1.intensity, 2.duration, 3.frequency, 4.autocorrelation,
%5.non-drought distribution
w = [1,4,1,2,2];
w = w/sum(w);
param.w = w;
param.tolerance = 0.02;


%% parameters used for experiment 3 "main_multi_drought.m"
% param.min_drought_duration  = 32; % minimum drought duration
% param.min_drought_intensity = -0.6; % minimum drought intensity
% param.ssi_time_scale        = 12; % aggregation time in months for SSI calculation 
% param.nmonth                = 60; % portion of the time series that is initially replaced in months
% param.nmonths_end_drought   = 12; % number of consecutive non-negative SSI values to end a drought
% 
% param.nyear_togenerate     = 100; % length of the generated scenarios in years
% param.H                    = param.nyear_togenerate*12; % % length of the generated scenarios in months
% 
% % define annealing schedule. 
% param.Nm = 1000; % number of swaps for each temperature and nmonth
% param.m  = 15 ; % number of temperature and nmonth reductions
% param.T  = 0.0001; % initial temperature
% param.decrease_rate = 0.8; %decrease rate of temperature and nmonth
% 
% %define weights 1.intensity, 2.duration, 3.frequency, 4.autocorrelation,
% %5.non-drought distribution
% w = [1,6,1,1,1];
% w = w/sum(w);
% param.w = w;
% param.tolerance = 0.005;
