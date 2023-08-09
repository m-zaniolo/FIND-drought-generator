function [Ointensity, Oduration, On_droughts, Oautocorr, Onon_drought_distrib] = ...
    calculate_objectives(obs, synth_flow, min_intensity, min_duration, time_scale, nmonths_end_drought,...
    target_int, target_pers, target_ndroguths, autocorr_hist, non_drought_distrib_hist)
% This function calculates the search objectives defined as follows:
% Ointensity: mean deviation between the intensity of each drought occurrence
%             in the synthetic scenario and the desired drought intensity
% Opersistence: mean deviation between the duration of each drought occurrence
%             in the synthetic scenario and the desired drought duration
% On_droughts: penalty term that activates when the number of droughts in the synthetic 
%              scenario is different than the desired number of droughts. 
% Oautocorr: mean deviation between the 1-year autocorrelation of
%            the historical time series and the synthetic one.
% Onon_drought_distrib: mean deviation between the 25th, 50th, and 75th
%           percentiles of non-drought periods in the historical and
%           synthetic time series.


% calculate drought index and statistics for synthetic flow time series
[duration_all, int_all, ndroguths_all, spi, drought_start_end] = drought_identification(...
    obs, synth_flow, min_intensity, min_duration, time_scale, nmonths_end_drought);  

%% calculate objectives
Ointensity  = mean(abs([int_all, mean(int_all)] - target_int));
Oduration = mean(abs([duration_all; mean(duration_all)] - target_pers));

if ndroguths_all < target_ndroguths % if number of droughts in synthetic time series is lower than target
On_droughts = 100*(target_ndroguths - ndroguths_all); % penalty proportional to the number of missing droughts 
elseif ndroguths_all > target_ndroguths % if number of droughts in synthetic time series is higher than target
  On_droughts = mean(duration_all); %penalty defined as average drought duration to guide search towards reducing drought periods
else
  On_droughts = 0; % no penalty if number of droughts in synthetic time series is equal to target
end

Oautocorr = mean(abs( autocorr_hist - correlogram(synth_flow, synth_flow, 12))); 

non_drought_periods = ones(length(spi),1);
% find non drought periods
for i = 1:size(drought_start_end, 1)
    non_drought_periods(drought_start_end(i, 1): drought_start_end(i, 2)) = 0;
end
% find 25th, 50th, and 75th percentile during non drought periods
non_drought          = spi(logical(non_drought_periods));
non_drought_distrib  = [prctile(non_drought, 25), median(non_drought), prctile(non_drought, 75) ];

Onon_drought_distrib = mean(abs(non_drought_distrib - non_drought_distrib_hist));

