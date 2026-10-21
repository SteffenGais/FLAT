function [flattened, baseline, slope] = flatten_spectrum(power, freqs, notch_hz)
%FLATTEN_SPECTRUM Flattening spectra using the Local Aperiodic Trend (FLAT).
% [flattened, baseline, slope] = flatten_spectrum(power, freqs)
% [flattened, baseline, slope] = flatten_spectrum(power, freqs, notch_hz)
% power: linear power/PSD, one row per spectrum and one column per frequency.
%        A row or column vector is accepted; its orientation is preserved.
% freqs: finite, nonnegative, strictly increasing frequencies in Hz.
% notch_hz: optional [low high] interval(s) to exclude, one per row.
%           Omitted or [] means no notch.
% flattened: power ./ baseline (dimensionless); baseline: input power units.
% Both are NaN outside 0.3–45 Hz and inside notch intervals. All-NaN rows
% remain NaN. Other analyzed power outside notches must be finite and positive.
% slope: one signed log-log baseline slope per spectrum over the full 2–24 Hz
% range, with equal log-frequency weighting. A 1/f^2 baseline gives -2.
% It is NaN if the full interval is unavailable or intersects a notch.
% The two-octave window and 1/32-octave grid are fixed method constants.
% Supply averaged spectra from a reliable passband; see docs/method.md.

if nargin < 3
    notch_hz = [];
end
original_size = size(power);
power = double(power);
if ~isvector(freqs)
    error('FLAT:Frequency', 'freqs must be a vector.');
end
freqs = double(freqs(:)');
if numel(freqs) < 2 || any(~isfinite(freqs)) || any(freqs < 0) || any(diff(freqs) <= 0)
    error('FLAT:Frequency', 'freqs must be finite, nonnegative, and strictly increasing.');
end
if isvector(power)
    power = power(:)';
end
if ~ismatrix(power) || isempty(power) || size(power, 2) ~= numel(freqs)
    error('FLAT:Shape', 'power must have one column per frequency.');
end
if ~isempty(notch_hz) && (size(notch_hz, 2) ~= 2 || any(~isfinite(notch_hz), 'all') ...
        || any(notch_hz(:, 1) <= 0) || any(notch_hz(:, 2) <= notch_hz(:, 1)) ...
        || any(notch_hz(2:end, 1) <= notch_hz(1:end-1, 2)))
    error('FLAT:Notch', 'Notch bands must be positive, ordered, disjoint [low high] pairs.');
end

% Work only in the fixed EEG/MEG range; never take the logarithm of DC.
selected = freqs >= 0.3 & freqs <= 45;
f = freqs(selected);
if numel(f) < 2
    error('FLAT:Range', 'At least two frequency bins are required within 0.3–45 Hz.');
end
x = log2(f);
grid = x(1):1/32:x(end)+1e-9;
if numel(grid) < 2
    error('FLAT:Range', 'The analyzed span must cover at least 1/32 octave.');
end
window = 65;                             % Two octaves: 64 intervals, 65 centered points.
half = floor(window / 2);
n_fit = min(half, numel(grid));           % Outer one octave at each edge.
offsets = -half:half;
[i, j] = find(triu(true(n_fit), 1));
edge_x = 0:n_fit-1;
masked = false(size(f));
for k = 1:size(notch_hz, 1)
    masked = masked | (f >= notch_hz(k, 1) & f <= notch_hz(k, 2));
end
if any(masked) && (masked(1) || masked(end))
    error('FLAT:NotchEdge', 'Each analyzed notch needs a frequency bin on both sides.');
end

flattened = NaN(size(power));
baseline = NaN(size(power));
slope = NaN(size(power, 1), 1);
for r = 1:size(power, 1)
    if all(isnan(power(r, :)))
        continue
    end
    p = power(r, selected);
    if any(~isfinite(p(~masked))) || any(p(~masked) <= 0)
        error('FLAT:Power', 'Analyzed power outside notches must be finite and positive.');
    end

    % Bridge excluded bins using unaffected neighbors in log-log coordinates.
    % Interpolated bins support the estimate but are never reported as data.
    logp = interp1(x(~masked), log10(p(~masked)), x, 'linear');
    % Roundoff can put the last grid point just beyond the measured endpoint.
    y = interp1(x, logp, min(grid, x(end)), 'linear');

    % Robust edge continuation prevents shrinking-window bias. Theil–Sen:
    % median pairwise slope, then median residual for the intercept.
    extensions = cell(1, 2);
    edges = {y(1:n_fit), y(end-n_fit+1:end)};
    queries = {-window:-1, n_fit:n_fit+window-1};
    for k = 1:2
        edge = edges{k};
        b = median((edge(j) - edge(i)) ./ (j' - i'));
        a = median(edge - b * edge_x);
        extensions{k} = a + b * queries{k};
    end
    y = [extensions{1}, y, extensions{2}];

    % First median-then-mean pass estimates the smooth local slope.
    first = movmean(movmedian(y, window), window);
    local_slope = gradient(first);

    % Remove the local slope BEFORE the median, so a positive peak cannot
    % shift the median's rank into a higher level of a sloping background.
    padded = [NaN(1, half), y, NaN(1, half)];
    neighborhoods = padded((1:numel(y))' + (0:window-1));
    detrended = neighborhoods - local_slope' .* offsets;
    trend = median(detrended, 2, 'omitnan')';
    trend = movmean(trend, window);
    trend = trend(window+1:window+numel(grid));

    % Continue the final segment if the last original bin falls just beyond
    % the regular grid; interp1's extrapolation avoids endpoint clamping.
    logbase = interp1(grid, trend, x, 'linear', 'extrap');
    b = 10 .^ logbase;
    b(masked) = NaN;
    baseline(r, selected) = b;
    flattened(r, selected) = p ./ b;

    % Equal log-frequency weighting, always over the full requested 2–24 Hz.
    overlaps_fit = ~isempty(notch_hz) && any(notch_hz(:, 1) <= 24 & notch_hz(:, 2) >= 2);
    if f(1) <= 2 && f(end) >= 24 && ~overlaps_fit
        fit_x = linspace(1, log2(24), 116);
        fit_y = interp1(x, logbase, fit_x, 'linear');
        centered = fit_x - mean(fit_x);
        slope(r) = sum(centered .* fit_y) / sum(centered .^ 2) / log10(2);
    end
end
flattened = reshape(flattened, original_size);
baseline = reshape(baseline, original_size);
end
