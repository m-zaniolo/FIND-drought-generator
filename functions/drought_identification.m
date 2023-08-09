function [persistence, intensity, frequency, index , drought_start_end] = drought_identification(data_ref, data_to_calc, min_drought_intensity, min_drought_duration, ssi_time_scale, nmonths_end_drought)

Ny     = length(data_ref)/12; 
year   = reshape( repmat(1:Ny, 12,1), Ny*12, 1); 
month  = repmat([1:12]', Ny, 1 ); 

method_name = 'SPI'; % the same method is used for both SSI and SPI computation  
% define your time scale here 
time_scale = ssi_time_scale; 
input_ref = [year, month, data_ref]; % prepare the input matrix 

Ny     = length(data_to_calc)/12; 
year   = reshape( repmat(1:Ny, 12,1), Ny*12, 1); 
month  = repmat([1:12]', Ny, 1 ); 
input_calc = [year, month, data_to_calc]; % prepare the input matrix 

sri = calc_drought_index(input_ref, time_scale, method_name, input_calc);
sri.index_value(sri.index_value<-4) = -4;
sri.index_value(sri.index_value>4)  = 4;

index = sri.index_value;
[persistence, intensity, frequency, drought_start_end ] = find_drought_periods(index, nmonths_end_drought, min_drought_duration, min_drought_intensity); 



