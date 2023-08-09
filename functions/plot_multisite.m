function plot_multisite(ssi_obs, ssi_obs2, ssi_synth, ssi_synth2, date, param)

hist_color   = [19, 108, 173]/256; % blue
synth_color  = [185, 76, 69]/256;  % red
obs2_color   = [1, 0.8398, 0];     % gold
synth2_color = [0 0.5 0];          % green

f = figure;
f.Position = [100 100 1250 400];
subplot(2,14,[1:7])
date_ssi = date(param.ssi_time_scale:end);
plot(date_ssi, ssi_obs ,Color=hist_color, LineWidth=1); hold on;
plot(date_ssi, ssi_obs2,Color=obs2_color, LineWidth=1); hold on;
plot(date_ssi, zeros(1, length(ssi_obs) ), 'k--');
legend('Site 1', 'Site 2');
title('Observed SSI');
ylim([-3, 3]);
set(gca, 'FontSize', 14); box on; grid on;

update_plot(1)

n_scenarios = size(ssi_synth, 2);
if n_scenarios>1
    b = uicontrol('Parent',f,'Style','slider','SliderStep', [1/(n_scenarios-1),1/(n_scenarios-1)], ...
        'Position',[81,54,23,320],'value',1, 'min',1, 'max',n_scenarios, 'Callback', @sliderCallback);

    % Add some text labels for the slider
    bgcolor = f.Color;
    bl1 = uicontrol('Parent',f,'Style','text','Position',[60,54,23,23],...
                    'String','1','BackgroundColor',bgcolor);
    bl2 = uicontrol('Parent',f,'Style','text','Position',[60,350,23,23],...
                    'String',num2str(n_scenarios),'BackgroundColor',bgcolor);
    bl3 = uicontrol('Parent',f,'Style','text','Position',[30,200,50,20],...
                    'String','scenario','BackgroundColor',bgcolor);

end

% Nested function to update the plot when the slider is moved
function update_plot(value)

synth  = ssi_synth(:,value);
synth2 = ssi_synth2(:,value);

subplot(2,14,[15:23]); cla();

date_synth_ssi = datetime(0,1,1)+calmonths(0:length(synth)-1);
plot(date_synth_ssi , synth ,Color=synth_color, LineWidth=1); hold on;
plot(date_synth_ssi , synth2,Color=synth2_color, LineWidth=1); hold on;
plot(date_synth_ssi , zeros(1, length(ssi_synth) ), 'k--');
legend('Site 1', 'Site 2');
title(['Synthetic SSI: Scenario ', num2str(value) ]);
ylim([-3, 3]);
set(gca, 'FontSize', 14); box on; grid on;

subplot(2,14,[11:14, 25:28]); cla();
visualize_deviation(ssi_obs, ssi_obs2, 0);
visualize_deviation(synth, synth2, 1);
corr_hist = corrcoef(ssi_obs,ssi_obs2);
corr_synth = corrcoef(synth,synth2);
legend(['observed data, \rho = ', num2str(corr_hist(1,2),2)] , ['synthetic data, \rho = ', num2str(corr_synth(1,2), 2)]);

function visualize_deviation(site1, site2, flag)

if flag == 0 % observed data
    color = hist_color;
    leg   = 'observed data';
else
    color = synth_color;
    leg   = 'synthetic data';
end


% Plot the scatterplot
scatter(site1, site2, 10, "filled", MarkerFaceColor=color); hold on;


xlim([-3, 3]); ylim([-3, 3]);
% Add axis labels and title
xlabel('Site 1');
ylabel('Site 2');
title('Multisite correlation');
set(gca, 'FontSize', 14); box on; grid on;
end
end

function sliderCallback(hObject, EventData)
    Value = round(get(hObject, 'Value'));
    set(hObject, 'Value', Value);
    update_plot(Value)

end
end
