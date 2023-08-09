function [processedData, dateVec_new] = data_preprocess(data_old, dateVec_old, time_scale)
% function processedData = data_preprocess(data_old)
% i dati in ingresso sono mensili e vengono processati per rimuovere il
% rumore, è inoltre possibile aggregarli in dati stagionali o comunque di
% più mesi 

% Input:
%            data_old - monthly time-series data
%         dateVec_old - date vector, with first column representing the
%                       year, and the second column for month
%          time_scale - scalar value for aggregating monthly
%                       time-series into sub-seasonal value
% Output:
%       processedData - processed time-series data 
%         dateVec_new - new date vector


no_years = unique(dateVec_old(:,1));
monthly_data_complete = nan(12, length(no_years));

no_missing_value = 12 - mod(length(data_old), 12);
if no_missing_value ~= 12
  warning(['Input time series is not dividable by 12. Missing month values',...
    ' will be filled with NaN.']);
end

for i = 1: length(no_years)
  idx1 = dateVec_old(:, 1) == no_years(i);
  idx2 = ismember(1:12, dateVec_old(idx1, 2));
  monthly_data_complete(idx2, i) = data_old(idx1);
end

if time_scale == 1
  processedData = monthly_data_complete(:);
  dateVec_new = dateVec_old;
  return
end

rolling_window = ones(time_scale, 1)/time_scale;

% Funzione FILTER (b,a,x): rimuove il rumore filtrando i dati secondo la funzione descritta
% da b e a. Qui è usato il filtro a media mobile, cioè a=1 e b=[1/m, 1/m, ...]

% notice here the NaN value will cause the processed Data
processedData = filter(rolling_window, 1, monthly_data_complete(:));
processedData(1: time_scale-1) = [];



dateVec_new = dateVec_old(time_scale:end, :);

end