"""Flatten the power spectra of a sleep recording per sleep stage.

    python sleep_example.py RECORDING SCORING [CHANNEL ...]

or in Python:

    from sleep_example import sleep_example
    result = sleep_example("night.vhdr", "night.txt", channels=["C3", "C4"])

RECORDING is a BrainVision (.vhdr), EDF, or other file MNE can read. SCORING
is a text file with one row per 30-s epoch, starting at the first sample of
the recording. Its first column holds the sleep stage: 0 wake, 1-4 NREM
stages 1-4 (3 and 4 are combined into N3), 5 REM. Other codes (e.g. 8 =
movement time) and further columns are ignored.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # flatten_spectrum.py
from flatten_spectrum import flatten_spectrum
from private.plots import plot_spectra, plot_waterfall
from private.recording import read_eeg
from private.spectra import log_cells

STAGES = ["W", "N1", "N2", "N3", "REM"]


def sleep_example(eeg_file, scoring_file, channels=None, frequency_range=(0.3, 45), notch_hz=None):
    """Average power spectra per sleep stage, flatten and plot them.

    One or two channels are shown as line plots, more as one waterfall figure
    per stage. channels, frequency_range and notch_hz work as in eeg_example.

    Returns a dict with freqs, power, baseline, flattened (power / baseline)
    and slope (2-24 Hz baseline slope) per stage and channel, and the
    channels, stages and number of epochs per stage.
    """
    raw = read_eeg(eeg_file, channels)
    fs = raw.info["sfreq"]
    epoch_samples = round(30 * fs)
    scores = np.loadtxt(scoring_file, ndmin=2)
    n_epochs = min(len(scores), raw.n_times // epoch_samples)

    # Stage of each epoch: 0 W, 1 N1, 2 N2, 3 N3, 4 REM; -1 = not analyzed.
    # Scoring codes 0 1 2 3 4 5 become W N1 N2 N3 N3 REM.
    stage_of_code = {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 4}
    stage = np.array([stage_of_code.get(code, -1) for code in scores[:n_epochs, 0]])
    # stage[scores[:n_epochs, 1] > 0] = -1      # also skip epochs marked in column 2 (MA, MT)
    if not np.any(stage >= 0):
        raise ValueError(f"No epochs with sleep stages 0-5 in {scoring_file}.")

    # Power spectrum of each epoch: one 30-s Hann window, power averaged within
    # 1/32-octave cells (full resolution at low, smooth spectra at high frequencies).
    spectra = {}
    for e in np.flatnonzero(stage >= 0):
        x = raw.get_data(start=e * epoch_samples, stop=(e + 1) * epoch_samples) * 1e6   # µV
        f, p = welch(x, fs, nperseg=epoch_samples, detrend="linear")
        keep = (f >= frequency_range[0]) & (f <= frequency_range[1])
        freqs, spectra[e] = log_cells(f[keep], p[:, keep])
        if (e + 1) % 100 == 0:
            print(f"Epoch {e + 1} of {n_epochs}")

    # Average the epochs of each stage in linear power, then flatten.
    epochs = np.array([np.sum(stage == s) for s in range(len(STAGES))])
    power = np.full((len(STAGES), len(raw.ch_names), freqs.size), np.nan)
    for s in np.flatnonzero(epochs):
        power[s] = np.mean([spectra[e] for e in np.flatnonzero(stage == s)], axis=0)
    power[np.any(power <= 0, axis=-1)] = np.nan     # flat channels have no spectrum
    flattened, baseline, slope = flatten_spectrum(power, freqs, notch_hz)

    name = Path(eeg_file).stem
    titles = [f"{stage_name} ({n} epochs)" for stage_name, n in zip(STAGES, epochs)]
    shown = np.flatnonzero(epochs)
    if len(raw.ch_names) <= 2:
        plot_spectra(freqs, power[shown], baseline[shown], flattened[shown], raw.ch_names,
                     [titles[s] for s in shown])
    else:
        for s in shown:
            plot_waterfall(freqs, power[s], flattened[s], raw.ch_names, f"{name}, {titles[s]}")
    return dict(freqs=freqs, power=power, baseline=baseline, flattened=flattened, slope=slope,
                channels=raw.ch_names, stages=STAGES, epochs=epochs)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    sleep_example(sys.argv[1], sys.argv[2], sys.argv[3:] or None)
    plt.show()
