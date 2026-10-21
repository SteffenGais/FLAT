function plot_spectra(freqs, power, baseline, flattened, labels, titles)
%PLOT_SPECTRA Line plots for one or two channels: power with its aperiodic
%   baseline (top) and the flattened spectrum, power / baseline (bottom).
%   power, baseline, flattened: channels x frequencies x conditions; one
%   column of plots per condition, titled by the string array titles.
n = size(power, 3);
figure('Color', 'w', 'Position', [100 100 200 + 260 * n, 600]);
tiledlayout(2, n, 'TileSpacing', 'compact');
colors = [0 0.447 0.741; 0.85 0.325 0.098];
top = gobjects(1, n);
bottom = gobjects(1, n);
for c = 1:n
    top(c) = nexttile(c);
    hold on
    for k = 1:size(power, 1)
        plot(freqs, power(k, :, c), 'Color', colors(k, :), 'LineWidth', 1.2, ...
            'DisplayName', labels{k});
        plot(freqs, baseline(k, :, c), '--', 'Color', colors(k, :), 'HandleVisibility', 'off');
    end
    set(gca, 'YScale', 'log');
    title(titles(c), 'Interpreter', 'none');

    bottom(c) = nexttile(n + c);
    hold on
    yline(1, ':', 'Color', [0.5 0.5 0.5]);
    for k = 1:size(power, 1)
        plot(freqs, flattened(k, :, c), 'Color', colors(k, :), 'LineWidth', 1.2);
    end
    xlabel('Frequency (Hz)');
end
linkaxes(top);
linkaxes(bottom);
set([top, bottom], 'XScale', 'log', 'XLim', freqs([1 end]), ...
    'XTick', [0.5 1 2 5 10 20 40], 'Box', 'on');
ylabel(top(1), 'Power (\muV^2/Hz)');
ylabel(bottom(1), 'Power / baseline');
plot(top(1), NaN, NaN, 'k--', 'DisplayName', 'Baseline');
legend(top(1), 'Location', 'southwest');
end
