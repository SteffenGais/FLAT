function plot_waterfall(freqs, power, flattened, labels, plot_title)
%PLOT_WATERFALL Waterfall plots of power (left) and the flattened spectrum,
%   power / baseline (right), for many channels. Channels with standard
%   10-05 names run from anterior (front) to posterior (back).
%   power, flattened: channels x frequencies.
order = 1:numel(labels);
positions = readtable(fullfile(fileparts(mfilename('fullpath')), 'standard_1005.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
[known, k] = ismember(labels, positions.name);
if all(known)
    [~, order] = sort(positions.y(k), 'descend');
end
keep = freqs <= 30;
figure('Color', 'w', 'Position', [50 100 1600 620]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
draw_waterfall(nexttile, freqs(keep), 10 * log10(power(order, keep)), labels(order), 'Power (dB)');
title('Power spectrum');
draw_waterfall(nexttile, freqs(keep), flattened(order, keep), labels(order), 'Power / baseline');
title('Flattened spectrum');
sgtitle(plot_title, 'Interpreter', 'none');
end


function draw_waterfall(ax, freqs, z, labels, z_label)
[F, C] = meshgrid(freqs, 1:numel(labels));
surf(ax, F, C, z, 'MeshStyle', 'row', 'LineWidth', 0.1);
colormap(ax, turbo);
view(ax, 9, 50);
box(ax, 'on');
grid(ax, 'on');
step = ceil(numel(labels) / 12);
set(ax, 'XLim', [0 30], 'YLim', [1 numel(labels)], ...
    'YTick', 1:step:numel(labels), 'YTickLabel', labels(1:step:end));
xlabel(ax, 'Frequency (Hz)');
zlabel(ax, z_label);
end
