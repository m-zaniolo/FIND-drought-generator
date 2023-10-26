function plot_objectives(obs, synthetic_flow, param)

n_scenarios = size(synthetic_flow, 2);
synth_color = [185, 76, 69]/256;  % red
hist_color = [19, 108, 173]/256;  % blue

f = figure;
f.Position = [100 100 250 1400];
% intensity
subplot(511); hold on;
plot([-1:n_scenarios+1], repmat(param.target_intensity, 1, n_scenarios+3), 'k--' )

% duration
subplot(512); hold on;
plot([-1:n_scenarios+1], repmat(param.target_duration, 1, n_scenarios+3), 'k--' )

% number of droughts
subplot(513); hold on;
plot([-1:n_scenarios+1], repmat(param.target_ndroughts, 1, n_scenarios+3), 'k--' )

% autocorrelation
subplot(514); hold on;
autocorr_hist = correlogram( obs, obs, 12 );
plot([0:12], autocorr_hist, LineWidth=2, Color=hist_color);
%plot(zeros(length(autocorr_hist),1), 'k');

% percentiles in historical non-drought periods
subplot(515); hold on;

[~, ~, ~, ssi_hist, drought_start_end] = drought_identification(obs, ...
        obs, param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale, param.nmonths_end_drought, param.distribution);

non_drought_periods = ones(length(ssi_hist), 1);
for i = 1:size(drought_start_end, 1)
    non_drought_periods(drought_start_end(i, 1): drought_start_end(i, 2)) = 0;
end
% find 25th, 50th, and 75th percentile during non drought periods
non_drought          = ssi_hist(logical(non_drought_periods));
non_drought_distrib_h  = [prctile(non_drought, 25), median(non_drought), prctile(non_drought, 75) ];

plot([-1:n_scenarios+1], repmat(non_drought_distrib_h(1), 1, n_scenarios+3), 'k--')
plot([-1:n_scenarios+1], repmat(non_drought_distrib_h(2), 1, n_scenarios+3), 'k--')
plot([-1:n_scenarios+1], repmat(non_drought_distrib_h(3), 1, n_scenarios+3), 'k--')


% plot historical drought properties
[duration_all, int_all, ndroguths_all, ssi_synth, drought_start_end] = drought_identification(obs, ...
    obs, param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale, param.nmonths_end_drought, param.distribution);

subplot(511);
scatter(zeros(size(int_all)), int_all, 'filled', MarkerFaceColor=hist_color, MarkerEdgeColor='none')

subplot(512);
scatter(zeros(size(duration_all)), duration_all, 'filled', MarkerFaceColor=hist_color, MarkerEdgeColor='none')

subplot(513);
bar(0, ndroguths_all+1, FaceColor=hist_color, EdgeColor="none")

subplot(515);
scatter(0, non_drought_distrib_h, 'filled', MarkerFaceColor=hist_color, MarkerEdgeColor='none')

% calculate and plot and objectives in all scenarios
for i = 1:n_scenarios
    [duration_all, int_all, ndroguths_all, ssi_synth, drought_start_end] = drought_identification(obs, ...
        synthetic_flow(:,i), param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale, param.nmonths_end_drought, param.distribution);

        subplot(511);
        scatter(i, int_all, 'filled', MarkerFaceColor=synth_color, MarkerEdgeColor='none')

        subplot(512);
        scatter(i, duration_all, 'filled', MarkerFaceColor=synth_color, MarkerEdgeColor='none')

        subplot(513);
        bar(i, ndroguths_all, FaceColor=synth_color, EdgeColor="none")

        subplot(514);
        autocorr_synth = correlogram( synthetic_flow(:,i), synthetic_flow(:,i), 12 );
        plot([0:12], autocorr_synth, Color=synth_color )

        subplot(515);
        non_drought_periods = ones(length(ssi_hist), 1);
        for j = 1:size(drought_start_end, 1)
        non_drought_periods(drought_start_end(j, 1): drought_start_end(j, 2)) = 0;
        end
        % find 25th, 50th, and 75th percentile during non drought periods
        non_drought          = ssi_synth(logical(non_drought_periods));
        non_drought_distrib  = [prctile(non_drought, 25), median(non_drought), prctile(non_drought, 75) ];

        scatter([i,i,i], non_drought_distrib, 'filled', MarkerFaceColor=synth_color, MarkerEdgeColor='none')

end

% adjust plot visual aspects

subplot(511);
set(gca, 'FontSize', 9); box on; grid on;
title('Drought intensity');
legend('Target intensity',  'Historical', 'Synthetic', Location='best');
xticks(0:n_scenarios)
xlim([-1, n_scenarios+1])
xticklabels(['H', strsplit(num2str(1:10))])
xlabel('Scenario')
ylabel('intensity [-]')
ylim([param.target_intensity + param.target_intensity*0.5 , param.target_intensity - param.target_intensity*0.5]);

subplot(512);
set(gca, 'FontSize', 9); box on; grid on;
title('Drought duration')
legend('Target duration',  'Historical', 'Synthetic', Location='best');
xticks(0:n_scenarios)
xlim([-1, n_scenarios+1])
xticklabels(['H', strsplit(num2str(1:10))])
xlabel('Scenario')
ylabel('duration [months]')
ylim([param.target_duration - param.target_duration*0.5 , param.target_duration + param.target_duration*0.5]);

subplot(513);
set(gca, 'FontSize', 9); box on; grid on;
title('Drought occurrences in 100y (Frequency)')
legend('Target frequency')
xticks(0:n_scenarios)
xlim([-1, n_scenarios+1])
xticklabels(['H', strsplit(num2str(1:10))])
xlabel('Scenario')
ylabel('number of droughts [-]')
ylim([ 0, param.target_ndroughts + 2]);

subplot(514);
plot([0:12], autocorr_hist, LineWidth=2, Color=hist_color); % so it's on top.
set(gca, 'FontSize', 9); box on; grid on;
title('Autocorrelogram')
xlim([0, 12])
ylabel('autocorrelation [-]')
legend('Historical autocorrelation')
xlabel('Lag')

subplot(515);
set(gca, 'FontSize', 9); box on; grid on;
xticks(0:n_scenarios)
xlim([-1, n_scenarios+1])
xticklabels(['H', strsplit(num2str(1:10))])
title('Percentiles in non drought periods')
legend('Target Percentiles', Location='best')
ylabel('SSI [-]')
xlabel('Scenario')
ylim([non_drought_distrib_h(1)*1.4, non_drought_distrib_h(3)*1.4])
