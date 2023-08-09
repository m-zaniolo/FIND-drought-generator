%% FIND drought generator

clc; clear;
addpath functions

rng(0); %reset random number generator for reproducibility

%% load data
obs  = load('Big Bend_1945_2022.txt');
date = datetime(1945,1,1)+calmonths(0:length(obs)-1)';

% Cut time series to retain only full years of data
obs  = obs(1:floor(length(obs)/12)*12 );
date = date(1:floor(length(obs)/12)*12 );
Ny   = length(date)/12;

%% define parameters for drought index calculation and synthetic time series optimization
param = generate_param_structure(); % this function contains important FIND parameters that can be modified by user
param.flag_visualize_optimization = 0; % set to 1 to visualize the ongoing progression of the search objective. Set to 0 otherwise.

%% define user-set drought parameters
% extract historical drought characteristics
[duration_obs, intensity_obs, n_droughts_obs, ssi_obs] = drought_identification(obs, obs, ...
    param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale, param.nmonths_end_drought);

% include in param the desired drought parameters
param.target_duration   = mean(duration_obs);
param.target_intensity  = mean(intensity_obs);
param.target_ndroughts  = round( n_droughts_obs/Ny *param.nyear_togenerate) ;
param.n_scenarios       = 10;

% some checks on data quality and selected characteristics 
check_drought_characteristics(obs, param)

%% experiment 1: run synthetic drought generator
[synthetic_flow, synthetic_ssi, duration_sim, intensity_sim] = run_synthetic_drought_generator(obs, param);

% plot results
plot_scenarios(date, obs, synthetic_flow, param);
plot_objectives(obs, synthetic_flow, param);

%% experiment 2: generate correlated streamflow scenarios for a second site

%load streamflow in site number 2
obs2  = load('Canby_1932_2022.txt');
date2 = datetime(1932,1,1)+calmonths(0:length(obs2)-1)';

% Cut time series to retain only full years of data
obs2  = obs2(1:floor(length(obs2)/12)*12);
date2 = date2(1:floor(length(obs2)/12)*12);

% intersect the 2 timeseries to obtain timespan when both sites
% measurements are available
[~, idx1, idx2] = intersect(date, date2);
date = date(idx1);
obs = obs(idx1);
date2 = date2(idx2);
obs2 = obs2(idx2);

% calculate SSI of site 2
[~, ~, ~, ssi_obs2] = drought_identification(obs2, obs2, ...
    param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale, param.nmonths_end_drought);
ssi_obs = ssi_obs( 1 + length(ssi_obs) - length(ssi_obs2):length(ssi_obs) );

% identify correlation between the 2 sites, defined as sum of squared
% difference between the 2. 
param.target_dispersion = sum((ssi_obs2 - ssi_obs ).^2); 

%run drought generator in multisite matching mode. 
[synthetic_flow2, synthetic_ssi2] = run_matching_drought(obs2, synthetic_ssi, param );
plot_multisite(ssi_obs, ssi_obs2, synthetic_ssi, synthetic_ssi2, date, param);

