function out = rolling_window_sum(vec, window_size)
% Computes the rolling window sum of a vector with a given window size
% Inputs:
%   - vec: the input vector
%   - window_size: the size of the window over which to compute the sum
% Output:
%   - out: the vector of rolling window sums

% Check for invalid window size
if window_size < 1 || window_size > length(vec)
    error('Invalid window size');
end

% Initialize the output vector with zeros
out = NaN*ones(size(vec));


% Compute the rolling window sum for the rest of the vector
for i = 1:length(vec)-window_size
    out(i) = sum(vec(i:i+window_size));
end

%out = out(1:end-window_size+1);

end
