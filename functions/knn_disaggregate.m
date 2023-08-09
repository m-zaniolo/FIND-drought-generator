function [disaggregated] = knn_disaggregate(cumulative, N, historical, idmonth)
% This function disaggregates a cumulative precipitation value into
% N monthly values using a k-NN approach with k=3 and historical monthly
% precipitation values.

% Determine the length of the historical precipitation vector
L = length(historical);

% Check that N is valid
if N < 1 || N > L
    error('Invalid value of N. N must be between 1 and the length of the historical precipitation vector.')
end

cum_historical = rolling_window_sum(historical, N);
cum_historical_r = reshape(cum_historical, 12, [])';

cum_historical_m = cum_historical_r(:, idmonth);

% Calculate the differences between the cumulative value and the historical
% values to find the nearest neighbors.
differences = abs(cumulative - cum_historical_m(1:sum(~isnan(cum_historical_m) ))  );

% Sort the historical values by their differences to the cumulative value
[sorted_differences, sorted_indices] = sort(differences);

% Use the k-NN approach to calculate the disaggregated values
k = 3; % Set the value of k to 3
neighbors = sorted_indices(1:k);

% extract one of the NNs with a probability proportional to their distance
inv_dists = 1 ./ (sorted_differences(1:k)+0.0001); % avoid NaN when distance is 0
inv_dists = inv_dists / sum(inv_dists);


idx = randsample(k, 1, true, inv_dists); 
selected_neighbor = neighbors(idx);

% adjust monthly values so their sum matches the cumulative value
nearest_indices_m = (selected_neighbor-1)*12+idmonth;
nearest_values = historical(nearest_indices_m:nearest_indices_m+N-1);
fractions = nearest_values./sum(nearest_values);

% fix if historical values are all null
if isnan(sum(fractions))
    fractions = (1/length(nearest_values))*ones(length(nearest_values), 1) ;
end

disaggregated = cumulative*fractions;

