function LMoment = pwm_Unbiased(data_to_calc)

[num_obs, num_samples] = size(data_to_calc);
data_to_calc_sort = sort(data_to_calc);

num_moments = 5;
[betas, L, R] = deal(nan(num_moments, num_samples));

for idx_col = 1: num_samples
  for r = 0: num_moments - 1
    %-- Calculate probability-weighted moments (PWM) values --
    idx_row = r+1;
    tmpVar = arrayfun( @nchoosek_modified, (0:num_obs-1)', r(ones(num_obs,1), :) )...
      .* data_to_calc_sort(:, idx_col);
    betas(idx_row, idx_col) = sum(tmpVar)/num_obs/nchoosek(num_obs-1, r);
    
    %-- Convert probability-weighted moments (PWM) to L-moments --
    weightSum = 0;
    for k = 0: r
      weight = (-1)^(r-k) * nchoosek(r, k) * nchoosek(r+k, k);
      weightSum = weightSum + weight * betas(k+1);
    end
    L(idx_row, idx_col) = weightSum;
  end
end

if (num_moments >= 2)
  R(2, :) = L(2, :)./L(1, :);
end

if (num_moments >= 3)
  R(3:end, :) = L(3:end, :)./L(2, :);
end

LMoment.lambda1 = L(1,:)';
LMoment.lambda2 = L(2,:)';
LMoment.lambda3 = L(3,:)';
LMoment.lambda4 = L(4,:)';
LMoment.lambda5 = L(5,:)';

LMoment.tau3 = R(3,:)';
LMoment.tau4 = R(4,:)';
LMoment.tau5 = R(5,:)';

end