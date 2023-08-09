function cdf_pearsonIII = calc_cdf_pearsonIII(param_pearsonIIIPdf, data_to_calc)

[num_obs, num_samples] = size(data_to_calc);
cdf_pearsonIII = nan(num_obs, num_samples);

if ( num_samples ~= length(param_pearsonIIIPdf.mu) )
  error('The number of samples does not match the number of parameters.')
end


SMALL = 1e-6;
for i = 1: num_samples
  mu = param_pearsonIIIPdf.mu(i);
  sigma = param_pearsonIIIPdf.sigma(i);
  gamma = param_pearsonIIIPdf.gamma(i);
  
  if (abs(gamma) <= SMALL)
    z_score = (data_to_calc(:, i) - mu) / sigma;
    cdf_pearsonIII(:, i) = normcdf(z_score);
  else
    ALPHA = 4/gamma^2;
    BETA = 0.5*sigma*abs(gamma);
    XI = mu - 2*sigma/gamma;
    
    if gamma > 0
      cdf_pearsonIII(:, i) = gamcdf((data_to_calc-XI)/BETA, ALPHA);
    else
      cdf_pearsonIII(:, i) = 1 - gamcdf((XI-data_to_calc)/BETA, ALPHA);
    end
    
  end
  
end

end