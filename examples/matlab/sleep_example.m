function result = sleep_example(eeg_file, scoring_file, options)
%SLEEP_EXAMPLE Flatten the power spectra of a sleep recording per stage.
%   SLEEP_EXAMPLE(EEG_FILE, SCORING_FILE) reads a BrainVision (.vhdr) or
%   EDF recording and a text file with one row per 30-s epoch, starting at
%   the first sample of the recording. The first column holds the sleep
%   stage: 0 wake, 1-4 NREM stages 1-4 (3 and 4 are combined into N3),
%   5 REM. Other codes (e.g. 8 = movement time) and further columns are
%   ignored. Power spectra are averaged per stage and flattened. One or two
%   channels are shown as line plots, more as one waterfall figure per stage.
%
%   Optional name-value arguments 'Channels', 'FrequencyRange' and 'Notch'
%   work as in EEG_EXAMPLE.
%
%   RESULT = SLEEP_EXAMPLE(...) returns the frequencies, power, baseline,
%   flattened spectrum and 2-24 Hz baseline slope per channel and
%   stage, and the number of epochs per stage.
%
%   Example:
%     sleep_example('night.vhdr', 'night.txt', 'Channels', {'C3', 'C4'})
arguments
    eeg_file {mustBeTextScalar}
    scoring_file {mustBeTextScalar}
    options.Channels = "all"
    options.FrequencyRange (1, 2) double = [0.3 45]
    options.Notch double = []
end
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));   % flatten_spectrum.m

rec = read_eeg(eeg_file, options.Channels);
epoch_samples = round(30 * rec.fs);
scores = readmatrix(scoring_file, 'FileType', 'text');
n_epochs = min(size(scores, 1), floor(rec.n_samples / epoch_samples));

% Stage of each epoch: 1 W, 2 N1, 3 N2, 4 N3, 5 REM; 0 = not analyzed.
% Scoring codes 0 1 2 3 4 5 become W N1 N2 N3 N3 REM.
names = ["W", "N1", "N2", "N3", "REM"];
stage_of_code = [1 2 3 4 4 5];
code = scores(1:n_epochs, 1);
stage = zeros(n_epochs, 1);
valid = ismember(code, 0:5);
stage(valid) = stage_of_code(code(valid) + 1);
% stage(scores(1:n_epochs, 2) > 0) = 0;      % also skip epochs marked in column 2 (MA, MT)
if ~any(stage)
    error('No epochs with sleep stages 0-5 in %s.', scoring_file);
end

% Power spectrum of each epoch: one 30-s Hann window, power averaged within
% 1/32-octave cells (full resolution at low, smooth spectra at high frequencies).
spectra = cell(n_epochs, 1);
for e = find(stage > 0)'
    x = rec.read((e - 1) * epoch_samples + 1, e * epoch_samples);
    [p, f] = welch_psd(x, rec.fs, 30);
    keep = f >= options.FrequencyRange(1) & f <= options.FrequencyRange(2);
    [freqs, spectra{e}] = log_cells(f(keep), p(:, keep));
    if mod(e, 100) == 0
        fprintf('Epoch %d of %d\n', e, n_epochs);
    end
end

% Average the epochs of each stage in linear power, then flatten.
n_stages = numel(names);
epochs = zeros(1, n_stages);
power = NaN(numel(rec.labels), numel(freqs), n_stages);
baseline = NaN(size(power));
flattened = NaN(size(power));
slope = NaN(numel(rec.labels), n_stages);
for s = 1:n_stages
    epochs(s) = sum(stage == s);
    if epochs(s) > 0
        p = mean(cat(3, spectra{stage == s}), 3);
        p(any(p <= 0, 2), :) = NaN;          % flat channels have no spectrum
        power(:, :, s) = p;
        [flattened(:, :, s), baseline(:, :, s), slope(:, s)] = ...
            flatten_spectrum(p, freqs, options.Notch);
    end
end

[~, name] = fileparts(eeg_file);
titles = compose("%s (%d epochs)", names', epochs');
shown = find(epochs > 0);
if numel(rec.labels) <= 2
    plot_spectra(freqs, power(:, :, shown), baseline(:, :, shown), ...
        flattened(:, :, shown), rec.labels, titles(shown));
else
    for s = shown
        plot_waterfall(freqs, power(:, :, s), flattened(:, :, s), rec.labels, ...
            sprintf('%s, %s', name, titles(s)));
    end
end

if nargout > 0
    result = struct('freqs', freqs, 'power', power, 'baseline', baseline, ...
        'flattened', flattened, 'slope', slope, 'channels', {rec.labels}, ...
        'stages', names, 'epochs', epochs);
end
end
