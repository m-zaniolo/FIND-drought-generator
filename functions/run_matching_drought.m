function [synth_flow, synth_ssi, obj] = run_matching_drought(obs, ssi_to_match_all, param)
% Generate synthetic time series of precipitation whose SSI
% has a specified dispersion wrt a given SSI
% % Inputs.
% obs:              historical flow time series in site 2
% ssi_to_match_all: matrix containing all scenarios of SSI generated for site 1
% param:            struct containing parameters for drought index calculation and synthetic time series optimization 
%
% Outputs.
% synth_flow:       generated streamflow time series for site 2
% synth_ssi:        SSI calculated for generated streamflow time series for site 2
% obj:              objective function value



%% Unpack parameters
min_intensity       = param.min_drought_intensity; % minimum drought duration
min_duration        = param.min_drought_duration; % minimum drought intensity
time_scale          = param.ssi_time_scale; % aggregation time in months for SSI calculation 
nyear_togenerate    = param.nyear_togenerate; % length of the generated scenarios in years
nmonths_end_drought = param.nmonths_end_drought; % number of consecutive non-negative SSI values to end a drought
nyear_inrecord      = length(obs)/12;  % length of the observed time series 


obs1_m = reshape(obs, 12, nyear_inrecord);

% search parameters
n_scenarios      = param.n_scenarios; 
Nm               = param.Nm; % number of swaps for each temperature
m                = param.m; % number of temperature and nmonth reductions
T                = param.T; % initial temperature
nmonth           = param.nmonth; % initial number of months consecutive months to replace
decrease_rate    = param.decrease_rate; % rate of temperature and nmonth decrease for every m

target_dispersion = param.target_dispersion;


% preallocate
O = zeros(n_scenarios,Nm*m-m);
simonth = zeros(n_scenarios,nyear_togenerate*12);

if param.flag_visualize_optimization
    figure;
end
%% loop through number of scenarios
for hh = 1:n_scenarios

    ssi_to_match = ssi_to_match_all(:,hh);
    duration_all = [];
    while isempty(duration_all) % check there is at least one drought in the initial random extraction
        Ny = nyear_togenerate; 
        os = [];
        acs = [];
        % for every month, fit data to distribution and extract random value
        for month = 1:12
          u = rand(1, Ny);
          LMoments = pwm_Unbiased(obs1_m(month,:)');
          param_pearsonIIIPdf = calc_param_pearsonIII(LMoments);
          accPdf = calc_cdf_pearsonIII(param_pearsonIIIPdf, obs1_m(month,:)');
    
          os(:, month)  = sort(obs1_m(month,:)');
          acs(:, month) = sort(accPdf);
          [~, un] = unique(acs(:, month));
          un_month{month} = un;
          xr(month, :) = interp1(acs(un, month),os(un, month),u,[],'extrap');
    
        end
    
        synth_flow = reshape(xr, [], 1);
        [duration_all] = drought_identification(obs, synth_flow, min_intensity, min_duration, time_scale, nmonths_end_drought, param.distribution); % check there is at least one drought in the initial random extraction
    end
    %% SIMULATED ANNEALING

    % preallocate matrices and set initialized values
    g               = repmat(synth_flow',[2 1]); % flow histogram
    Objective       = nan;
    k               = 0;
    should_continue = 1; 

    while should_continue % check if termination criteria are met
        
        % algorithm starts at high temperature and then decreases so that you are less likely to accept a nonimproving swap
        k = k+1;
        T = decrease_rate.*T;

        % algorithm starts by swapping a large number of months and progressively 
        % decreases to refine scenarios at a shorter time scale
        nmonth = max(1, round( decrease_rate.*nmonth) );
        
        a = 1; % lower bound random sampling
        b1 = nyear_togenerate*12 - nmonth; % upper bound random sampling

        % cumulate historical data to the timescale of the iteration 
        cum_prec = rolling_window_sum(obs, nmonth);
        cum_prec_m = reshape(cum_prec, 12, nyear_inrecord);
        cum_prec_m = cum_prec_m(:, ~isnan( sum(cum_prec_m) ) );

        % fit cumulative data to distribution 
        os = [];
        acs = [];

        for month = 1:12

          cum_prec = cum_prec_m( ~isnan( cum_prec_m(month,:)'));
          LMoments = pwm_Unbiased(cum_prec);
          param_pearsonIIIPdf = calc_param_pearsonIII(LMoments);
          accPdf = calc_cdf_pearsonIII(param_pearsonIIIPdf, cum_prec);
    
          os(:, month)  = sort(cum_prec);
          acs(:, month) = sort(accPdf);
          [~, un] = unique(acs(:, month));
          un_month{month} = un;
    
        end


        for i = 2:Nm %number of iterations for each temperature
            if should_continue

                %% calculate statistics of old and swaped time series
    
                if i == 2 % in first iteration, compute SSI
                  [~, ~, ~, ssi_prev] = drought_identification(obs, g(1,:)', min_intensity, min_duration, time_scale, nmonths_end_drought, param.distribution);  % returns intensity and persistence of all the detected droughts
    
                elseif as == 1 % if swap in previous iteration is accepted
                    ssi_prev   = ssi_swap;
                end
    
    
                %% swap new value into synthetic time series
                
                % extract a random element from the synthetic time series to be
                % the starting point for the swap
                r = round(a + (b1-a).*rand(1,1));
                
                idmonth = rem(r, 12);
                if idmonth == 0
                    idmonth = 12;
                end
    
                % extract random value from the relative cumulated monthly distribution
                u = rand(1);
                cum_value = interp1(acs(un_month{idmonth}, idmonth),os(un_month{idmonth}, idmonth),u,[],'extrap');
                if cum_value < 0
                    cum_value = 0;
                end
                % disaggregate cumulated value into monthly
                disagg = knn_disaggregate(cum_value, nmonth, obs, idmonth);
    
                % assign this value to swapped time series
                g(2,r:r+nmonth-1) = disagg ; 
    
                %% calculate statistics after swapping
                  [~, ~, ~, ssi_swap] = drought_identification(obs, g(2,:)', min_intensity, min_duration, time_scale, nmonths_end_drought, param.distribution); 
                  
    
                %% calculate value of aggregate objective function for swapped and unchanged cases
    
                dispersion_prev = sum((ssi_prev - ssi_to_match).^2); % dispersion is defined as sum of squared differences between synthetic SSI and SSI to match
                dispersion_swap  = sum((ssi_swap - ssi_to_match).^2); 
    
                OFold  = abs(target_dispersion - dispersion_prev); % dispersion should match the target dispersion
                OFswap = abs(target_dispersion - dispersion_swap); 
    
                pmov = min( 1, exp(((OFold-OFswap)/OFold)/T)); % probability of swapping the two values
    
                %do the swapping
                if OFswap<=OFold % swaps that reduce the value of the objective function are kept
                    %keep swap
                    g(1,:) = g(2,:); % swap values
                    as = 1;
                else % swaps that increase the value are  kept according to the pmov value
                    if rand<pmov % generate a random number. random number uniformly distributed between 0 and 1.
                        % we select one such number and compare it with pmov.
                        % keep swap
                        g(1,:) = g(2,:); % swap values
                        as = 1;
                    else
                        % no swap
                        g(2,:)    = g(1,:);
                        as = 0;
                    end
                end
    
                Objective = [Objective; OFold];
                if param.flag_visualize_optimization
                    plot(Objective); 
                    drawnow;
                end
                if OFold < param.tolerance % termination criterion on objective tolerance
                    should_continue = 0;
                end
            end
        end

       if k >= m %interrupt when the max number of iterations is reached
           should_continue = 0;
       elseif OFold<30 % or, interrupt when objective function is good enough. 
           should_continue = 0;
       end

    end %%

%% log results of terminated simulation

    O(hh,1:length(Objective)-1 ) = Objective(2:end);
    ext_sri = g(1,:)';
    simonth(hh,:) = ext_sri';


end

%% generate outputs
obj   = O(:,end);
synth_flow  = simonth';

for hh = 1:n_scenarios
    [~, ~, ~, synth_ssi(:,hh)] = drought_identification(obs, synth_flow(:,hh), min_intensity, min_duration, time_scale, nmonths_end_drought, param.distribution);
end
