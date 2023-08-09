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

