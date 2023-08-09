clc; clear;
addpath functions

rng(0); %reset random number generator for reproducibility

%% load data
obs  = load('Big Bend_1945_2022.txt');
date = datetime(1945,1,1)+calmonths(0:length(obs)-1)'; 

% Cut time series to eliminate any spare months that don't belong to a full
% year of data
obs  = obs(1:floor(length(obs)/12)*12);
date = date(1:floor(length(obs)/12)*12);
Ny = length(date)/12;

%% define parameters for drought index calculation and synthetic time series optimization
param = generate_param_structure();
param.flag_visualize_optimization = 0; % set to 1 to visualize the ongoing progression of the search objective. Set to 0 otherwise.


%% experiment 3: generate a set of scenarios with varyinf drought characteristics of intensity and duration
% extract historical drought characteristics
[duration_obs, intensity_obs, n_droughts_obs, ssi_obs] = drought_identification(obs, obs, param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale, param.nmonths_end_drought);

ccD = [0.75, 1, 1.25, 1.5, 1.75]; %drought duration multiplier. 
ccI = [0.75, 1, 1.25, 1.5, 1.75]; %drought intensity multiplier. 
param.n_scenarios       = 3;      % number of scenarios to generate
k = 0;
synth_flow_all = [];
synth_ssi_all = [];

% loop through defined drought characteristics
for i = 1:length(ccD)
    for j = 1:length(ccI)
        k = k+1;
        % include in param the desired drought parameters
        param.target_duration   = mean(duration_obs)*ccD(i);
        param.target_intensity  = mean(intensity_obs)*ccI(j);
        param.target_ndroughts  = n_droughts_obs ;
        
        % some checks on data quality and selected characteristics 
        check_drought_characteristics(obs, param)
        
        %% run synthetic drought generator
        [synthetic_flow, synthetic_ssi, duration_sim, intensity_sim, obj] = run_synthetic_drought_generator(obs, param);
        synth_flow_all = [synth_flow_all, synthetic_flow];
        synth_ssi_all  = [synth_ssi_all, synthetic_ssi];
    end
end

%% plot matrix of drought duration and intensity and two selected trajectories

plot_multidrought(synth_flow_all, obs, param, duration_obs, intensity_obs, ccD, ccI);
