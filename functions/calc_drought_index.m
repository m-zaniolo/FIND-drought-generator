function results_summary = calc_drought_index(data_ref, time_scale, indicator_name, varargin)
%RESULTS_SUMMARY   control function for drought index toolbox
%
%   CALC_DROUGHT_INDEX(data_ref, time_scale, indicator_name, varargin)
%   calculate the drought index given the input data and specified method
%
%
%  Example: 
%  
%  Author: Yu Li (yu.li@polimi.it)
%  Date: 30th March, 2016

% ---- error check ----

m = size(data_ref);
if ( m~=3 )
  error('The input data must have three columns, i.e., [year, month, value]')
end

% ---- configure the computation route ----
data_to_calc = data_ref; 

if nargin > 3 
    data_to_calc = varargin{1};
end
if nargin > 4
    distribution = varargin{2};
else 
    distribution = 'PearsonIII';
end
    
    
if ~strcmpi(distribution, 'gamma')
    if ~strcmpi(distribution, 'PearsonIII')
        fprintf('Undefined distribution. Using Pearson III instead')
        distribution = 'PearsonIII';
    end
end

if strcmpi(indicator_name, 'SPI')
  fooHandle = @calc_SPI;
  UseMethod = distribution;

end

%fprintf(1, ['Calculating the %s index now. Estimating the parameters', ...
%' assuming ''%s'' distribution! \n'], upper(indicator_name), upper(UseMethod));
    
% ---- pre-process the dataset based on the specified time scale ----

[data_ref_processed, idx_date_ref] = data_preprocess(...
  data_ref(:,3), data_ref(:,1:2), time_scale);
[data_to_calc_processed, idx_date_calc] = data_preprocess(...
  data_to_calc(:,3), data_to_calc(:,1:2), time_scale);
  
% ---- start calculation ----

Mnths = unique(idx_date_calc(:,2));
index_values = zeros(size(idx_date_calc, 1), 1);

for i = 1: length(Mnths)
  eachMnth = Mnths(i);
  idx1 = idx_date_ref(:, 2) == eachMnth;
  idx2 = idx_date_calc(:,2) == eachMnth;
  
  index_values(idx2) = fooHandle( data_ref_processed(idx1), UseMethod, data_to_calc_processed(idx2) );
end

% ---- summarize the results ----

results_summary.time_scale = time_scale;
results_summary.date_index = idx_date_calc;
results_summary.index_name.name = upper(indicator_name);
results_summary.index_name.method = upper(UseMethod);
results_summary.index_value = index_values;


end