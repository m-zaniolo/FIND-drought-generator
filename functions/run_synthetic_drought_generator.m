function [synth_flow, synth_ssi, duration_sim, intensity_sim, ndroughts, obj] = run_synthetic_drought_generator(obs, param)
% Generate synthetic time series with specified drought characteristics
% Inputs.
% obs:             historical flow time series
% param:           struct containing parameters for drought index calculation and synthetic time series optimization 
%
% Outputs.
% synth_flow:      generated streamflow time series
% synth_ssi:       SSI calculated for generated streamflow time series
% duration_sim:    duration of each drought in the generated streamflow time series 
% intensity_sim:   intensity of each drought in the generated streamflow time series 
% ndroughts:       number of droughts in the generated streamflow time series 
% obj:             objective function value

%% Unpack parameters
% parameters related to drought calculation and historical time series
min_intensity       = param.min_drought_intensity; % minimum drought duration
min_duration        = param.min_drought_duration; % minimum drought intensity
time_scale          = param.ssi_time_scale; % aggregation time in months for SSI calculation 
nyear_togenerate    = param.nyear_togenerate; % length of the generated scenarios in years
nmonths_end_drought = param.nmonths_end_drought; % number of consecutive non-negative SSI values to end a drought
nyear_inrecord      = length(obs)/12;  % length of the observed time series 


% search parameters
n_scenarios       = param.n_scenarios; % number of generated scenarios
Nm                = param.Nm; % number of swaps for each temperature
m                 = param.m; % number of temperature and nmonth reductions
T                 = param.T; % initial temperature
nmonth            = param.nmonth; % initial number of months consecutive months to replace
decrease_rate     = param.decrease_rate; % rate of temperature and nmonth decrease for every m

% target drought parameters
target_intensity = param.target_intensity; % desired drought intensity
target_duration  = param.target_duration; % desired drought duration
target_ndroughts = param.target_ndroughts; % desired number of droughts in the synthetic time series

w = param.w;


%% calculate historical observed statistics of interest

% observed monthly autocorrelation
autocorr_hist = correlogram( obs, obs, 12);

% Calculate the Standardized Streamflow Index (SSI) and the drought occurrence time of the input time series
[~, ~, ~, ssi_obs , drought_start_end_obs] = drought_identification(obs, obs, min_intensity, min_duration, time_scale, nmonths_end_drought, param.distribution);  % compute SRI from streamflow data
obs1_m      = reshape(obs, 12, nyear_inrecord);

% Find non-drought periods and calculate observed 25th, 50th, and 75th percentile
non_drought_periods_hist = ones(length(ssi_obs), 1);
for i = 1:size(drought_start_end_obs, 1)
    non_drought_periods_hist(drought_start_end_obs(i, 1): drought_start_end_obs(i, 2)) = 0;
end

non_drought_hist          = ssi_obs(logical(non_drought_periods_hist));
non_drought_distrib_hist  = [prctile(non_drought_hist, 25), median(non_drought_hist), prctile(non_drought_hist, 75) ];


% preallocate
O = zeros(n_scenarios,Nm*m-m);
simonth = zeros(n_scenarios,nyear_togenerate*12);


%% loop through number of scenarios
for hh = 1:n_scenarios

    duration_all = [];
    while isempty(duration_all) % check there is at least one drought in the initial random extraction

        Ny = nyear_togenerate; 
        os = [];
        acs = [];
        
        % for every month, fit monthly data to distribution and extract random value
        for month = 1:12
          u = rand(1, Ny);
          accPdf = fit_distribution(obs1_m(month,:)', param.distribution); % choose 'gamma' or 'Pearson III' distributions

          os(:, month)  = sort(obs1_m(month,:)');
          acs(:, month) = sort(accPdf);
          [~, un] = unique(acs(:, month));
          un_month{month} = un;
          xr(month, :) = interp1(acs(un, month),os(un, month),u,[],'extrap');

          

    
        end

        % obtain inizialized synthetic time series
        synth_flow = reshape(xr, [], 1); 

        % K-S test of synthetic data:
        isFit = kstest2(sort(obs) ,sort(synth_flow));
        if param.flag_visualize_optimization
            if hh == 1
                if isFit
                    fprintf('The K-S test fails to reject the null hypothesis at a default 5 percent significance level. The probability distribution selected to approximate the data is appropriate')
                else
                    fprintf('The K-S test rejects the null hypothesis at a default 5 percent significance level')
                end
            end
        end


        [duration_all] = drought_identification(obs, synth_flow, min_intensity, min_duration, time_scale, nmonths_end_drought, param.distribution); 
    end
    %% Optimization portion

    % preallocate matrices and set initialized values
    Objective       = nan;
    k               = 0;
    should_continue = 1; 

    % matrix g contains the previous (first row) and swapped (second row)
    % synthetic time series and is updated at every iteration. Both rows
    % are initialized here with the synth_flow generated above.
    g               = repmat(synth_flow',[2 1]); 

    
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
        cum_flow = rolling_window_sum(obs, nmonth);
        cum_flow_m = reshape(cum_flow, 12, nyear_inrecord);
        cum_flow_m = cum_flow_m(:, ~isnan( sum(cum_flow_m) ) );

        % fit cumulative data to distribution 
        os = [];
        acs = [];
        for month = 1:12

          cum_flow = cum_flow_m( ~isnan( cum_flow_m(month,:)'));
          LMoments = pwm_Unbiased(cum_flow);
          param_pearsonIIIPdf = calc_param_pearsonIII(LMoments);
          accPdf = calc_cdf_pearsonIII(param_pearsonIIIPdf, cum_flow);
    
          os(:, month)  = sort(cum_flow);
          acs(:, month) = sort(accPdf);
          [~, un] = unique(acs(:, month));
          un_month{month} = un;
    
        end


        for i = 2:Nm %number of iterations for each temperature
            if should_continue

                %% calculate objectives of previous time series
    
                if i == 2 % in first iteration, compute intensity, persistence, frequency pre- and post-swap
                  [Ointensity_prev, Oduration_prev, Ondroughts_prev, Oautocorr_prev, Onon_drought_distrib_prev] = calculate_objectives(obs, ...
                      g(1,:)', min_intensity, min_duration, time_scale, nmonths_end_drought, target_intensity, target_duration, target_ndroughts, ...
                      autocorr_hist, non_drought_distrib_hist, param.distribution);
    
                elseif as == 1 % if swap in previous iteration is accepted
                    Ointensity_prev           = Ointensity_swap;
                    Oduration_prev            = Oduration_swap;
                    Ondroughts_prev           = Ondroughts_swap;
                    Oautocorr_prev            = Oautocorr_swap;
                    Onon_drought_distrib_prev = Onon_drought_distrib_swap;
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
    
    
                %% calculate objectives of swapped time series
                  [Ointensity_swap, Oduration_swap, Ondroughts_swap, Oautocorr_swap, Onon_drought_distrib_swap] = calculate_objectives(obs, ...
                      g(2,:)', min_intensity, min_duration, time_scale, nmonths_end_drought, target_intensity, target_duration, target_ndroughts, ...
                      autocorr_hist, non_drought_distrib_hist, param.distribution);
    
    
                %% calculate partial and aggregate objective function for previous and swapped cases
                %a rescaling factor is introduced for the duration
                %objective so all objectives have comparable order of
                %magnitude
                OFduration_prev  = Oduration_prev/100;  
                OFduration_swap  = Oduration_swap/100;
    
                % objective aggregation
                OFold  = w(1)*Ointensity_prev + w(2)*OFduration_prev + w(3)*Oautocorr_prev + w(4)*Ondroughts_prev  + w(5)*Onon_drought_distrib_prev;
                OFswap = w(1)*Ointensity_swap + w(2)*OFduration_swap + w(3)*Oautocorr_swap + w(4)*Ondroughts_swap  + w(5)*Onon_drought_distrib_swap;
    
    
                % calculate probability of swapping the two values as a
                % function of temperature
                pmov = min( 1, exp(((OFold-OFswap)/OFold)/T)); 
    
                % do the swapping
                if OFswap<=OFold % swaps that reduce the value of the objective function are kept
                    %keep swap
                    g(1,:) = g(2,:); % swap values
                    as = 1;
                    %OFold = OFswap;
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
                    if OFold<10
                        ylim([0, 10])
                        if OFold<1
                            ylim([0, 1])
                            if OFold<0.2
                               ylim([0, 0.2])
                            end
                        end
                    end
                   drawnow;
                end 
                if OFold < param.tolerance % termination criterion on objective tolerance
                    should_continue = 0;
                end

            end
        end

       % termination criterion for max NFE
       if k >= m %interrupt when the max number of iterations is reached
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
    [duration, int, ndroughts(hh), synth_ssi(:,hh)] = drought_identification(obs, synth_flow(:,hh), min_intensity, min_duration, time_scale, nmonths_end_drought, param.distribution);
    duration_sim(hh)  = mean(duration);
    intensity_sim(hh) = mean(int);
end
