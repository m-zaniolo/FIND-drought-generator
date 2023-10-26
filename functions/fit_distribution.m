function accPdf = fit_distribution(data, useMethod)
  

if strcmpi(useMethod, 'PearsonIII')
    LMoments = pwm_Unbiased(data);
    param_pearsonIIIPdf = calc_param_pearsonIII(LMoments);
    accPdf = calc_cdf_pearsonIII(param_pearsonIIIPdf, data);
elseif strcmpi(useMethod, 'gamma')
    accPdf = calc_cdf_gamma(data);
else
    fprintf('Undefined distribution. Using PearsonIII instead')
    LMoments = pwm_Unbiased(obs1_m(month,:)');
    param_pearsonIIIPdf = calc_param_pearsonIII(LMoments);
    accPdf = calc_cdf_pearsonIII(param_pearsonIIIPdf, data);

end

end

function accPdf = calc_cdf_gamma(data)

logData = nan(size(data));
idxZero = data==0; % find the zero precipitation estimates to avoid INF after log operation

mu_ref = mean(data(~idxZero));
logData(~idxZero) = log(data(~idxZero));
U_val = log(mu_ref) - nanmean(logData);

shapeParameter = ( 1 + sqrt(1 + 4*U_val/3) )/(4*U_val);
scaleParameter = mu_ref / shapeParameter;

% -- Calculate the cumulative denstiy function --
accPdf = gamcdf(data, shapeParameter, scaleParameter);
end
