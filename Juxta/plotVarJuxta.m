function plotVarJuxta(allCJuxta, field)
% plotVarJuxta - Plots histogram and subtype distribution for a given field.
%
% Inputs:
%   allCJuxta - Struct: Contains data fields for Juxta analysis.
%   field     - String: Field name to plot.
%
% The function:
%   - Extracts and processes the field values.
%   - Plots a histogram.
%   - Categorizes and plots subtype distributions.

% Extract data from structure
si = getfield(allCJuxta, field);
identity = allCJuxta.identity;

% Process specific fields
if strcmp(field, 'pbAwakeChange') || strcmp(field, 'pbSleepChange')
    si = abs(si);
end

% Define histogram edges
edges = linspace(min(si), max(si), 20);

% Maximize figure window
set(gcf, 'WindowState', 'maximized');

% Histogram Plot
figure;
histogram(si, edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.5);
title(field);
set(gca, 'Color', 'none', 'Box', 'off', 'TickDir', 'out', 'FontSize', 15);
ylimVals = ylim;
yticks(0:floor(ylimVals(2)/3):ylimVals(2));
hold on;

% Subtype Labeling
labels = {'pv', 'som', 'pvsom', 'nothing'}; % Categories
colors = {'b', 'r', 'm', 'k'}; % Corresponding Colors
yOffset = 1;

% Plot individual subtype points
for a = 1:length(labels)
    subtypeValues = si(identity == labels{a});
    for i = 1:length(subtypeValues)
        plot(subtypeValues(i), -3 * yOffset * a, '.', 'Color', colors{a}, 'MarkerSize', 50);
        hold on;

        % Count duplicate occurrences and mark them
        count = sum(subtypeValues == subtypeValues(i));
        if count > 1
            plot(subtypeValues(i), -2 * yOffset * a - 0.1, '.', 'Color', colors{a}, 'MarkerSize', 50);
        end
    end
end

ylim([-16 0]);
set(gca, 'Color', 'none', 'Box', 'off', 'YColor', 'none', 'TickDir', 'out');

end
