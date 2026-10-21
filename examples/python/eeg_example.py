"""Flatten the power spectra of an EEG recording and plot them.

Without a file, one minute of 64-channel resting EEG with eyes closed is
downloaded (PhysioNet EEG Motor Movement/Imagery Dataset, subject 98).

    python eeg_example.py
    python eeg_example.py recording.vhdr
    python eeg_example.py recording.vhdr O1 O2

or in Python:

    from eeg_example import eeg_example
    result = eeg_example("recording.vhdr", channels=["O1", "O2"])

Data: Schalk et al. (2004), IEEE Trans Biomed Eng 51:1034-1043, and
Goldberger et al. (2000), Circulation 101:e215-e220 (license ODC-BY 1.0).
"""

import sys
from pathlib import Path
from urllib.request import urlretrieve

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import welch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # flatten_spectrum.py
from flatten_spectrum import flatten_spectrum
from private.plots import plot_spectra, plot_waterfall
from private.recording import read_eeg
from private.spectra import log_cells

DATA = Path(__file__).resolve().parents[1] / "data"
URL = "https://physionet-open.s3.amazonaws.com/eegmmidb/1.0.0/S098/S098R02.edf"


def eeg_example(eeg_file=None, channels=None, frequency_range=(0.3, 45), notch_hz=None):
    """Flatten the average power spectrum of each channel and plot it.

    eeg_file: BrainVision (.vhdr), EDF, or any other file MNE can read;
        None downloads the example recording.
    channels: channel names (default: all). One or two channels are shown as
        line plots with their aperiodic baseline, more channels as waterfalls.
    frequency_range: (low, high) in Hz that the recording measures reliably,
        i.e. away from filter edges (default: the full FLAT range).
    notch_hz: (low, high) bands to exclude, one row per band (default: none).

    Returns a dict with freqs, power, baseline, flattened (power / baseline),
    slope (2-24 Hz baseline slope) and channels.
    """
    if eeg_file is None:
        eeg_file = DATA / "S098R02.edf"
        if not eeg_file.exists():
            DATA.mkdir(parents=True, exist_ok=True)
            print(f"Downloading {eeg_file.name}")
            part = eeg_file.with_name(eeg_file.name + ".part")
            urlretrieve(URL, part)
            part.rename(eeg_file)

    # Average power spectrum over the whole recording, read in blocks of 30 s so
    # that long recordings fit into memory. Hann windows (50% overlap within a
    # block), as long as the recording allows while averaging about 30 of them
    # (4-30 s).
    raw = read_eeg(eeg_file, channels)
    fs = raw.info["sfreq"]
    seconds = min(30, max(4, raw.n_times / fs / 15))
    block = min(round(30 * fs), raw.n_times)
    n_blocks = raw.n_times // block
    power = 0
    for b in range(n_blocks):
        x = raw.get_data(start=b * block, stop=(b + 1) * block) * 1e6   # µV
        freqs, p = welch(x, fs, nperseg=round(seconds * fs), detrend="linear")
        power = power + p / n_blocks

    # Keep the reliable frequency range, average power within 1/32-octave cells
    # (full resolution at low, smooth spectra at high frequencies) and flatten.
    keep = (freqs >= frequency_range[0]) & (freqs <= frequency_range[1])
    freqs, power = log_cells(freqs[keep], power[:, keep])
    power[np.any(power <= 0, axis=-1)] = np.nan     # flat channels have no spectrum
    flattened, baseline, slope = flatten_spectrum(power, freqs, notch_hz)

    name = Path(eeg_file).stem
    if len(raw.ch_names) <= 2:
        plot_spectra(freqs, power[None], baseline[None], flattened[None], raw.ch_names, [name])
    else:
        plot_waterfall(freqs, power, flattened, raw.ch_names,
                       f"{name}, {len(raw.ch_names)} channels")
    return dict(freqs=freqs, power=power, baseline=baseline, flattened=flattened,
                slope=slope, channels=raw.ch_names)


if __name__ == "__main__":
    eeg_example(sys.argv[1] if len(sys.argv) > 1 else None, sys.argv[2:] or None)
    plt.show()
