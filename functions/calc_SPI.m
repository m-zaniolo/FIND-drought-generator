function index_SPI = calc_SPI(data_ref, useMethod, varargin)
%CALC_SPI  calculate the SPI index
%
%   CALC_SPI(data_ref) calculate the standard precipitation index (SPI)
%   given the reference precipitation estiamtes data_ref. The precipitation
%   estimates are fitted with Gamma distribution or Pearson Type III
%   distribution and translated into normal one, with zero mean. Positive
%   SPI value indicates a wet period, while a negative value indicates a
%   dry period.
%
%   CALC_SPI(data_ref, data_to_calc) calculate the SPI for prediction by
%   using data_ref for estimating the Gamma distribution, and data_to_calc
%   for computing the SPI value
%
%
%  Author: Yu Li (yu.li@polimi.it)
%  Date: 30th March, 2016
%
%  Reference:
%
%  Naresh Kumar, M., Murthy, C.S., Sesha Sai, M.V.R., Roy, P.S., 2009. On
%  the use of Standardized Precipitation Index (SPI) for drought intensity
%  assessment. Meteorological Applications 16, 381–389.

% -- Prepare the data --
data_to_calc = data_ref;

if (nargin > 2)
  data_to_calc = varargin{1};
end

if ~isvector(data_ref) || ~isvector(data_to_calc)
  error('All input data must be vector.');
end

no_obs = length(data_ref);

% -- Fit the data with Gamma or Pearson-III distribution --

if strcmpi(useMethod, 'gamma')
  accPdf = calc_cdf_gamma(data_ref, data_to_calc);
  
elseif strcmpi(useMethod, 'PearsonIII')  
  LMoments = pwm_Unbiased(data_ref);
  param_pearsonIIIPdf = calc_param_pearsonIII(LMoments);
  accPdf = calc_cdf_pearsonIII(param_pearsonIIIPdf, data_to_calc);
else
  error('Method undefined!')
end


% Since the Gamma function is undefined for x = 0 and a precipitation
% distribution may contain zeros, the cumulative probability needs to be
% adjusted
p_zero = sum(data_to_calc == 0)/no_obs; % probability of zero precipitation
accPdf_adjust = p_zero + (1-p_zero).*accPdf;

% -- Translate the gamma distribution into standard normal distribution --
index_SPI = norminv(accPdf_adjust, 0, 1);

idxInf = isinf(index_SPI);
signVal = sign(index_SPI(idxInf));

% Force the Inf/-Inf value to be 3/-3, indicating an extreme wet/dry
% situation
index_SPI(idxInf) = 3*signVal;


end



%% ========================================================================
% =========================================================================
% ========== Subfunctions =================================================

function accPdf = calc_cdf_gamma(data_ref, data_to_calc)

logData = nan(size(data_ref));
idxZero = data_ref==0; % find the zero precipitation estimates to avoid INF after log operation

mu_ref = mean(data_ref(~idxZero));
logData(~idxZero) = log(data_ref(~idxZero));
U_val = log(mu_ref) - nanmean(logData);

shapeParameter = ( 1 + sqrt(1 + 4*U_val/3) )/(4*U_val);
scaleParameter = mu_ref / shapeParameter;

% -- Calculate the cumulative denstiy function --
accPdf = gamcdf(data_to_calc, shapeParameter, scaleParameter);
end

function param_pearsonIIIPdf = calc_param_pearsonIII(LMoments)

SMALL = 1e-06;

C1  =  0.2906 ;
C2  =  0.1882 ;
C3  =  0.0442 ;
D1  =  0.36067;
D2  = -0.59567;
D3  =  0.25361;
D4  = -2.78861;
D5  =  2.56096;
D6  = -0.77045;
PI3 =  3 * pi ;

rootPi = sqrt(pi);

L1 = LMoments.lambda1;
L2 = LMoments.lambda2;
T3 = abs(LMoments.tau3);

param_pearsonIIIPdf = struct('mu',[], 'sigma',[], 'gamma',[]);
[param_pearsonIIIPdf.mu, param_pearsonIIIPdf.sigma, param_pearsonIIIPdf.gamma] = ...
  deal(zeros(size(L1)));

[T, ALPHA] = deal(nan(size(T3)));
idx = T3 >= 1/3;
T(idx) = 1-T3;
ALPHA(idx) = T(idx) .* (D1 + T(idx) .* (D2 + T(idx) .* D3))./...
  (1 + T(idx) .* (D4 + T(idx) .* (D5 + T(idx) .* D6)));

T(~idx) = PI3 * T3.^2;
ALPHA(~idx) = (1 + C1 .* T(~idx))./...
  (T(~idx) .* (1 + T(~idx) .* (C2 + T(~idx) .* C3)));

RTALPH = sqrt(ALPHA);
BETA = rootPi .* L2 .* exp(gammaln(ALPHA) - gammaln(ALPHA + 0.5));

param_pearsonIIIPdf.mu = L1;
param_pearsonIIIPdf.sigma = BETA .* RTALPH;
param_pearsonIIIPdf.gamma = 2/RTALPH;

idx = T3 <= SMALL;
param_pearsonIIIPdf.mu(idx) = L1;
param_pearsonIIIPdf.sigma(idx) = L2.* rootPi;
param_pearsonIIIPdf.gamma(idx) = 0;

idx = LMoments.tau3 < 0;
param_pearsonIIIPdf.gamma(idx) = -param_pearsonIIIPdf.gamma(idx);


end


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