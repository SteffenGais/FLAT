# FLAT for Python
[Overview](README.md) · [MATLAB](README_matlab.md) · [Method](docs/method.md)

## Installation

FLAT needs Python 3.10 or newer and NumPy. The examples also use SciPy,
Matplotlib and MNE.

With pip:

```sh
pip install flatten-spectrum               # FLAT only
pip install "flatten-spectrum[examples]"   # FLAT and the packages for the examples
```

With conda, create an environment and install with pip inside it:

```sh
conda create -n flat python=3.12 pip
conda activate flat
pip install "flatten-spectrum[examples]"
```

To install from a copy of this repository, run `pip install ".[examples]"` in
its main folder. You can also copy `flatten_spectrum.py` into your project; it
only needs NumPy.

## Usage

```python
from flatten_spectrum import flatten_spectrum

flattened, baseline, slope = flatten_spectrum(power, freqs)
flattened, baseline, slope = flatten_spectrum(power, freqs, notch_hz=[(47.5, 52.5)])
```

Inputs:

- `power`: power spectral density in linear units (e.g. µV²/Hz), with
  frequencies along the last axis. Leading axes (channels, conditions, ...)
  are kept. Preferably from long windows with 1/32-octave cells (see below).
- `freqs`: frequencies in Hz, increasing.
- `notch_hz` (optional): frequency bands to exclude, as `(low, high)` pairs.
  Default: none.

Outputs:

- `flattened`: the flattened spectrum, `power / baseline`. It is 1 at the
  background level and shows periodic activity as peaks, together with the
  noise of the estimate.
- `baseline`: the aperiodic background, in the units of `power`.
- `slope`: the slope of the baseline from 2 to 24 Hz in log-log coordinates,
  one per spectrum (−2 for a 1/f² background). The aperiodic exponent of
  specparam (FOOOF, fixed mode) corresponds to −slope.

Frequencies outside 0.3–45 Hz and inside notch bands are NaN in `flattened`
and `baseline`. Pass only the frequencies your recording measures reliably,
away from filter edges (see [method](https://github.com/SteffenGais/FLAT/blob/main/docs/method.md#reliable-frequency-range)):

```python
keep = (freqs >= 1) & (freqs <= 40)
flattened, baseline, slope = flatten_spectrum(power[..., keep], freqs[keep])
```

## Recommended spectra

Compute the spectrum from windows as long as the recording allows while
still averaging about 30 of them, e.g. one 30-s Hann window per sleep epoch,
or 4-s windows with 50% overlap for one minute of data. Average the
periodograms in linear power, then average power within 1/32-octave cells:
low frequencies keep their full resolution, high frequencies become smooth
(see [method](https://github.com/SteffenGais/FLAT/blob/main/docs/method.md#input-spectra)).

```python
import numpy as np
from scipy.signal import welch

freqs, power = welch(data, fs, nperseg=int(10 * fs))   # data: channels x samples, 10-s windows
cells = np.round(32 * np.log2(freqs[1:])).astype(int)  # without the DC bin
_, first, counts = np.unique(cells, return_index=True, return_counts=True)
freqs = np.add.reduceat(freqs[1:], first) / counts
power = np.add.reduceat(power[..., 1:], first, axis=-1) / counts
```

## Examples

The examples in `examples/python` read any recording MNE can read, e.g.
BrainVision (`.vhdr`) or EDF. Downloads are stored in `examples/data`.

```sh
cd examples/python
python eeg_example.py                           # 1 min of 64-channel resting EEG (1.3 MB download)
python eeg_example.py my_recording.vhdr O1 O2   # your own recording, two channels
python sleep_edf_example.py                     # one night of the Sleep-EDF database (50 MB download)
python sleep_example.py night.vhdr night.txt    # your own sleep recording and scoring
```

One or two channels are shown as line plots with their baseline, more
channels as waterfall plots of power and flattened power.

`sleep_example.py` needs a text file with one row per 30-s epoch, starting at
the first sample of the recording. Its first column holds the sleep stage:
0 wake, 1–4 NREM stages 1–4, 5 REM. Other codes (e.g. 8 for movement time) are
skipped, and stages 3 and 4 are combined into N3.

In Python or Jupyter, the examples are functions with the optional arguments
`channels`, `frequency_range` and `notch_hz`. They return their results:

```python
from sleep_example import sleep_example

result = sleep_example("night.vhdr", "night.txt", channels=["C3", "C4"])
```
