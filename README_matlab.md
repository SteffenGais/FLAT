# FLAT for MATLAB
[Overview](README.md) · [Python](README_python.md) · [Method](docs/method.md)

`flatten_spectrum.m` needs MATLAB R2018b or newer, the examples R2020b or
newer. No toolboxes are required.

## Setup

Add the FLAT folder to the MATLAB path:

```matlab
addpath('C:\path\to\FLAT')
```

## Usage

```matlab
[flattened, baseline, slope] = flatten_spectrum(power, freqs);
[flattened, baseline, slope] = flatten_spectrum(power, freqs, notch_hz);
```

Inputs:

- `power`: power spectral density in linear units (e.g. µV²/Hz), one row per
  spectrum (e.g. channel) and one column per frequency, preferably from long
  windows with 1/32-octave cells (see below).
- `freqs`: frequencies in Hz, increasing.
- `notch_hz` (optional): frequency bands to exclude, one `[low high]` row per
  band. Default: none.

Outputs:

- `flattened`: the flattened spectrum, `power ./ baseline`. It is 1 at the
  background level and shows periodic activity as peaks, together with the
  noise of the estimate.
- `baseline`: the aperiodic background, in the units of `power`.
- `slope`: the slope of the baseline from 2 to 24 Hz in log-log coordinates,
  one per row (−2 for a 1/f² background). The aperiodic exponent of
  specparam (FOOOF, fixed mode) corresponds to −slope.

Frequencies outside 0.3–45 Hz and inside notch bands are NaN in `flattened`
and `baseline`. Pass only the frequencies your recording measures reliably,
away from filter edges (see [method](docs/method.md#reliable-frequency-range)):

```matlab
keep = freqs >= 1 & freqs <= 40;
[flattened, baseline, slope] = flatten_spectrum(power(:, keep), freqs(keep));
```

## Recommended spectra

Compute the spectrum from windows as long as the recording allows while
still averaging about 30 of them, e.g. one 30-s Hann window per sleep epoch,
or 4-s windows with 50% overlap for one minute of data. Average the
periodograms in linear power, then average power within 1/32-octave cells:
low frequencies keep their full resolution, high frequencies become smooth
(see [method](docs/method.md#input-spectra)).

```matlab
% power: channels x frequencies, freqs: 1 x frequencies (transpose pwelch output)
f = freqs(freqs > 0);
p = power(:, freqs > 0);
[~, ~, cells] = unique(round(32 * log2(f)));
members = sparse(1:numel(cells), cells, 1);    % bins x cells
counts = full(sum(members, 1));
freqs = (f * members) ./ counts;
power = full(p * members) ./ counts;
```

The examples do this without toolboxes (using `welch_psd.m` and `log_cells.m` in
`examples/matlab/private`).

## Examples

The examples in `examples/matlab` read BrainVision (`.vhdr`) and EDF files.
Downloads are stored in `examples/data`.

```matlab
cd examples/matlab
eeg_example                                 % 1 min of 64-channel resting EEG (1.3 MB download)
eeg_example('', 'Channels', {'O1', 'O2'})   % two channels as line plots
eeg_example('my_recording.vhdr')            % your own recording
sleep_edf_example                           % one night of the Sleep-EDF database (50 MB download)
sleep_example('night.vhdr', 'night.txt')    % your own sleep recording and scoring
```

One or two channels are shown as line plots with their baseline, more
channels as waterfall plots of power and flattened power.

`sleep_example` needs a text file with one row per 30-s epoch, starting at the
first sample of the recording. Its first column holds the sleep stage: 0 wake,
1–4 NREM stages 1–4, 5 REM. Other codes (e.g. 8 for movement time) are
skipped, and stages 3 and 4 are combined into N3.

`eeg_example` and `sleep_example` accept the optional arguments `'Channels'`
(names or numbers), `'FrequencyRange'` (`[low high]` in Hz, default
`[0.3 45]`) and `'Notch'`, and return their results:

```matlab
result = sleep_example('night.vhdr', 'night.txt', 'Channels', {'C3', 'C4'});
```
