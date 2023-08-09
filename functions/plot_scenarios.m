function plot_scenarios(date, obs, synthetic_all, param )
% Generates a GUI to plot historical and synthetic climate data
% Inputs:
%   - date: vector of datetime values representing the time period of the data
%   - obs: vector of observed climate data
%   - synthetic_all: matrix of synthetic climate data, where each column
%     represents a different scenario
%   - param: structure containing model parameters
%
% Outputs:
%   - None (generates a plot)

% Define colors for plotting
hist_color = [19, 108, 173]/256;  % blue
synth_color = [185, 76, 69]/256;  % red

% Create a figure
f = figure;

% Determine the number of scenarios
n_scenarios = size(synthetic_all, 2);

% Plot historical SSI and droughts
subplot(2,14,[1:7])


[~, ~, n_droughts, ssi, drought_start_end] = drought_identification(obs, obs, param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale, param.nmonths_end_drought );
ssi = [zeros(param.ssi_time_scale-1,1); ssi]; %SSI timeseries is 11 month shorter than streamflow index because of aggregation
date_ssi = date;
area(date_ssi, ssi ,'EdgeColor','none', FaceColor=[175, 171, 171]/265); hold on;
for i = 1: n_droughts
    % Highlight the drought periods in red.
    a = drought_start_end(i,1)+param.ssi_time_scale-1;
    b = drought_start_end(i,2)+param.ssi_time_scale-1;
    area(date_ssi(a:b), ssi(a:b), 'FaceColor', 'r', 'EdgeColor','none' )
end
xlim([date(1), date(end)]);
ax1 = gca;
ax1.YLabel.String = 'SSI [-]';
ax1.YLim=[-3.5, 3.5];
% add streamflow timeseries
C = colororder;
yyaxis right
area(date_ssi, obs/1000, 'EdgeColor',C(5,:) )
yticks([0 5 10])
ylim([0, 25])
ylabel('flow [th cfs]                    ')

% Add title, legend, adjust figure appearence
title('Observed SSI')
legend('Normal and wet conditions', 'Droughts'); %, Location='eastoutside');

set(gca, 'FontSize', 14); box on; grid on;
f.Position = [100 100 1250 400];

% Call update_figure to initialize the plot with results from first
% scenario
update_figure(1)

% Add a slider if there are multiple scenarios
if n_scenarios>1
    b = uicontrol('Parent',f,'Style','slider','SliderStep', [1/(n_scenarios-1),1/(n_scenarios-1)], ...
        'Position',[81,54,23,320],'value',1, 'min',1, 'max',n_scenarios, 'Callback', @sliderCallback);

    % Add some text labels for the slider
    bgcolor = f.Color;
    bl1 = uicontrol('Parent',f,'Style','text','Position',[60,54,23,23],...
                    'String','1','BackgroundColor',bgcolor);
    bl2 = uicontrol('Parent',f,'Style','text','Position',[60,350,23,23],...
                    'String',num2str(n_scenarios),'BackgroundColor',bgcolor);
    bl3 = uicontrol('Parent',f,'Style','text','Position',[30,200,50,20],'fontsize',12,...
                    'String','scenario','BackgroundColor',bgcolor);

end

% Nested function to update the plot when the slider is moved
function update_figure(value)
    % If multiple scenarios are provided, select the scenario based on the
    % slider value. Otherwise, use the only scenario available.
    if n_scenarios > 1
        synthetic = synthetic_all(:,value);
    else
        synthetic = synthetic_all;
    end

    % plot synthetic drought scenario
    subplot(2,14,[15:23]); cla();
    yyaxis left; cla();
    
    

    % Identify droughts in the synthetic scenario
    [~, ~, frequency, ssi, drought_start_end] = drought_identification(obs, synthetic, param.min_drought_intensity, param.min_drought_duration, param.ssi_time_scale, param.nmonths_end_drought );
    ssi = [zeros(param.ssi_time_scale-1,1); ssi]; %SSI timeseries is 11 month shorter than streamflow index because of aggregation
    
    % Plot the SSI values for the synthetic scenario
    date_synth_ssi = datetime(0,1,1)+calmonths(0:length(ssi)-1);
    area(date_synth_ssi, ssi ,'EdgeColor','none', FaceColor=[175, 171, 171]/265); hold on;
    for i = 1: frequency
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
    title(['Synthetic SSI: Scenario ', num2str(value) ])
    %legend('Normal and wet conditions', 'Droughts');
    set(gca, 'FontSize', 14); box on; grid on;
    
    xlabel('Time')

    % add streamflow timeseries
    C = colororder;
    yyaxis right
    area(date_synth_ssi, synthetic/1000, 'EdgeColor',C(5,:), 'FaceColor',C(5,:) )
    yticks([0 5 10])
    ylim([0, 25])
    ylabel('flow [th cfs]                    ')
    
    ax2 = gca;
    ax2.YColor = C(5,:);



    % Compare cyclostationary monthly values of streamflow in historical and
    % synthetic scenario
    subplot(2,14,[11:14, 25:28]); cla();
    obs2 = obs/1000;
    % compute monthly statistics for historical data
    monthly = reshape(obs2, 12, [])' ;
    m = prctile(monthly, 25);
    M = prctile(monthly, 75);
    avg_hist = median(monthly);

    % plot historical data statistics
    x = 1:12;
    hold on;
    x2 = [x, fliplr(x)];
    inBetween = [m, fliplr(M)];
    fill(x2, inBetween, hist_color, 'FaceAlpha',0.5, 'LineStyle', 'none');

    % compute monthly statistics for synthetic data
    monthly = reshape(synthetic/1000, 12, [])' ;
    m = prctile(monthly, 25);
    M = prctile(monthly, 75);
    avg_synth = median(monthly);

    % plot synthetic data statistics
    inBetween = [m, fliplr(M)];
    fill(x2, inBetween, synth_color, 'FaceAlpha',0.5, 'LineStyle', 'none');
    plot(avg_hist, 'LineWidth',2, 'Color',hist_color);
    plot(avg_synth, 'LineWidth',2, 'Color',synth_color);

    % add title, legend, and axes labels
    legend('25th to 75th percentile historical', '25th to 75th percentile synthetic', 'Median historical', 'Median synthetic', Location='northeast')
    title(['Monthly streamflow, scenario: ', num2str(value)])
    ylabel('Flow [thousand cfs]')
    xticks([1:12])
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'});
    xlim([1,12])
    set(gca, 'FontSize', 14); box on; grid on;

end

function sliderCallback(hObject, EventData)
    Value = round(get(hObject, 'Value'));
    set(hObject, 'Value', Value);
    update_figure(Value)

end
end
