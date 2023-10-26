function plot_multidrought(synth_flow_all, obs, param, duration_obs, intensity_obs, ccD, ccI)

f = figure;

% setup bottom axis
ax_top = axes(); % axis to appear at top
set(gca, 'FontSize', 12);
hold(ax_top);
ax_top.XAxisLocation = 'top';
ax_top.YAxisLocation = "right";
ax_top.Color = 'none';
xlabel(ax_top, 'Duration multiplier [-]')
ylabel(ax_top, 'Intensity multiplier [-]')


ax = axes();
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = "left";
set(gca, 'FontSize', 12);
hold(ax);
xlabel(ax, 'Duration [months]'); 
ylabel(ax, 'intensity [-]');
% setup top axis



k = 0; 
c = [];
target_ = [];
for j = 1:5 
    for i = 1:5
        k = k+1;
        % define color scheme
        ii = 0.2:0.16:1;
        jj = 0.2:0.16:1;
        c(k,:) = rgb2hsv([1,1,1].*[ii(i), jj(j), 0.5]); 
        c(k,:) = [1,1,1].*[ii(i), 0, jj(j)];
        
        % target drought characteristics
        target(k,:)   = [mean(duration_obs)*ccD(i), mean(intensity_obs)*ccI(j)];
        target_(i,j) = k ;
    end
    
end

% define colormap
h = colormap(summer(5));
map = [];
weights = 0.5:0.15:1.15;
for i = 1:5
    map = [map; min(ones(size(h)), h*weights(i)) ];
end
c = map; 

% plot grid
target_ = [target_-1, zeros(5, 1); zeros(1, 6)];
X = [mean(duration_obs)*ccD - 0.125*mean(duration_obs), mean(duration_obs)*ccD(length(ccD)) + 0.125*mean(duration_obs) ]; 
Y = [mean(intensity_obs)*ccI  - 0.125*mean(intensity_obs), mean(intensity_obs)*ccI(length(ccI)) + 0.125*mean(intensity_obs) ]; 

a = pcolor(X, Y, target_); hold on ;

cc = colormap(gca, c);
a.FaceAlpha = 0.9;
for k = 1:25
    for kk = 1:param.n_scenarios
        [duration, intensity] = drought_identification(obs, synth_flow_all(:,(k-1)*param.n_scenarios + kk), param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale, param.nmonths_end_drought, param.distribution );
        scatter(duration, intensity, 60, cc(k,:), 'filled', 'MarkerEdgeColor', 'k'); %c(k,4)
    end
end

%plot target as white dot
TD = repmat(mean(duration_obs)*ccD, 1, 5);
TI = reshape( repmat(mean(intensity_obs)*ccI, 5, 1), 1, []);

scatter(TD,TI, 10, 'white', 'filled'); 

%set axis limit
half_step = (ccD(2) - ccD(1))/2;
ax_top.XLim=[min(ccD) - half_step, max(ccD) + half_step];
ax_top.XTick = ccD;

ax_top.YLim=[min(ccI) - 0.125, 1.75 + 0.125];
ax_top.YTick = ccI;
ax_top.YTickLabel = flip(ccI);

ax.XLim = [min(X), max(X)];
ax.YLim = [min(Y), max(Y)];

%% trajectory plots

%
f2 = figure;
f2.Position = [500 150 800 400];
subplot(211)
id_short_intense = length(ccI)*param.n_scenarios + 1; % index of drought scenario with shortest duration and highest intensity
prec = synth_flow_all(:, id_short_intense); 

[~, ~, n_droughts, ssi, drought_start_end] = drought_identification( obs,prec, param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale, param.nmonths_end_drought, param.distribution );
    ssi = [zeros(param.ssi_time_scale-1,1); ssi]; %SSI timeseries is ssi_time_scale month shorter than streamflow index because of aggregation
    
    % Plot the SSI values for the synthetic scenario
    date_synth_ssi = datetime(0,1,1)+calmonths(0:length(ssi)-1);
    area(date_synth_ssi, ssi ,'EdgeColor','none', FaceColor=[175, 171, 171]/265); hold on;
    for i = 1: n_droughts
        % Highlight the drought periods in red.
    a = drought_start_end(i,1)+param.ssi_time_scale-1;
    b = drought_start_end(i,2)+param.ssi_time_scale-1;
        area(date_synth_ssi(a:b), ssi(a:b), 'FaceColor','r', 'EdgeColor','none' )
    end
    ax = gca;
    ax.YColor = 'k';
    ax.YLabel.String = 'SSI [-]';
    ax.YLim=[-3.5, 3.5];
    % Set the title, legend, and figure appearence
    title('Short intense droughts')
    %legend('Normal and wet conditions', 'Droughts');
    set(gca, 'FontSize', 14); box on; grid on;
    
    xlabel('Time')

    % add streamflow timeseries
    C = colororder;
    yyaxis right
    area(date_synth_ssi, prec/1000, 'EdgeColor',C(5,:), 'FaceColor',C(5,:) )
    yticks([0 5 10])
    ylim([0, 25])
    ylabel('flow [th cfs]                    ')
        ax2 = gca;
    ax2.YColor = C(5,:);


subplot(212)
id_long_mild = length(ccI)*(length(ccD)-1)*param.n_scenarios + 1; % index of drought scenario with longest duration and mildest intensity
prec = synth_flow_all(:, id_long_mild);
[~, ~, n_droughts, ssi, drought_start_end] = drought_identification( obs,prec, param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale,param.nmonths_end_drought , param.distribution );
    ssi = [zeros(param.ssi_time_scale-1,1); ssi]; %SSI timeseries is 11 month shorter than streamflow index because of aggregation
    
    % Plot the SSI values for the synthetic scenario
    date_synth_ssi = datetime(0,1,1)+calmonths(0:length(ssi)-1);
    area(date_synth_ssi, ssi ,'EdgeColor','none', FaceColor=[175, 171, 171]/265); hold on;
    for i = 1: n_droughts
        % Highlight the drought periods in red.
    a = drought_start_end(i,1)+param.ssi_time_scale-1;
    b = drought_start_end(i,2)+param.ssi_time_scale-1;
        area(date_synth_ssi(a:b), ssi(a:b), 'FaceColor','r', 'EdgeColor','none' )
    end
    ax = gca;
    ax.YColor = 'k';
    ax.YLabel.String = 'SSI [-]';
    ax.YLim=[-3.5, 3.5];
    % Set the title, legend, and figure appearence
    title('Long mild droughts')
    %legend('Normal and wet conditions', 'Droughts');
    set(gca, 'FontSize', 14); box on; grid on;
    
    xlabel('Time')

    % add streamflow timeseries
    C = colororder;
    yyaxis right
    area(date_synth_ssi, prec/1000, 'EdgeColor',C(5,:), 'FaceColor',C(5,:) )
    yticks([0 5 10])
    ylim([0, 25])
    ylabel('flow [th cfs]                    ')
        ax2 = gca;
    ax2.YColor = C(5,:);



