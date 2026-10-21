function result = eeg_example(eeg_file, options)
%EEG_EXAMPLE Flatten the power spectra of an EEG recording and plot them.
%   EEG_EXAMPLE downloads one minute of 64-channel resting EEG with eyes
%   closed (PhysioNet EEG Motor Movement/Imagery Dataset, subject 98) and
%   plots its power spectra before and after flattening.
%
%   EEG_EXAMPLE(EEG_FILE) uses a BrainVision (.vhdr) or EDF (.edf) file.
%
%   Optional name-value arguments:
%     'Channels'        channel names or numbers (default: all). One or two
%                       channels are shown as line plots with their
%                       aperiodic baseline, more channels as waterfalls.
%     'FrequencyRange'  [low high] in Hz that the recording measures
%                       reliably, i.e. away from filter edges
%                       (default [0.3 45], the full FLAT range).
%     'Notch'           [low high] bands to exclude, one row per band
%                       (default: none).
%
%   RESULT = EEG_EXAMPLE(...) returns the frequencies, power, baseline,
%   flattened spectrum and 2-24 Hz baseline slope per channel.
%
%   Example:
%     eeg_example('', 'Channels', {'O1', 'O2'})
%
%   Data: Schalk et al. (2004), IEEE Trans Biomed Eng 51:1034-1043, and
%   Goldberger et al. (2000), Circulation 101:e215-e220 (license ODC-BY 1.0).
arguments
    eeg_file {mustBeTextScalar} = ""
    options.Channels = "all"
    options.FrequencyRange (1, 2) double = [0.3 45]
    options.Notch double = []
end
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));   % flatten_spectrum.m

if eeg_file == ""
    url = 'https://physionet-open.s3.amazonaws.com/eegmmidb/1.0.0/S098/S098R02.edf';
    eeg_file = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'S098R02.edf');
    if ~isfile(eeg_file)
        if ~isfolder(fileparts(eeg_file))
            mkdir(fileparts(eeg_file));
        end
        fprintf('Downloading S098R02.edf\n');
        websave([eeg_file '.part'], url);
        movefile([eeg_file '.part'], eeg_file);
    end
end

% Average power spectrum over the whole recording, read in blocks of 30 s so
% that long recordings fit into memory. Hann windows (50% overlap within a
% block), as long as the recording allows while averaging about 30 of them
% (4-30 s).
rec = read_eeg(eeg_file, options.Channels);
seconds = min(30, max(4, rec.n_samples / rec.fs / 15));
block = min(round(30 * rec.fs), rec.n_samples);
n_blocks = floor(rec.n_samples / block);
power = 0;
for b = 1:n_blocks
    [p, freqs] = welch_psd(rec.read((b - 1) * block + 1, b * block), rec.fs, seconds);
    power = power + p / n_blocks;
end

% Keep the reliable frequency range, average power within 1/32-octave cells
% (full resolution at low, smooth spectra at high frequencies) and flatten.
keep = freqs >= options.FrequencyRange(1) & freqs <= options.FrequencyRange(2);
[freqs, power] = log_cells(freqs(keep), power(:, keep));
power(any(power <= 0, 2), :) = NaN;          % flat channels have no spectrum
[flattened, baseline, slope] = flatten_spectrum(power, freqs, options.Notch);

[~, name] = fileparts(eeg_file);
if numel(rec.labels) <= 2
    plot_spectra(freqs, power, baseline, flattened, rec.labels, string(name));
else
    plot_waterfall(freqs, power, flattened, rec.labels, ...
        sprintf('%s, %d channels', name, numel(rec.labels)));
end

if nargout > 0
    result = struct('freqs', freqs, 'power', power, 'baseline', baseline, ...
        'flattened', flattened, 'slope', slope, 'channels', {rec.labels});
end
end
