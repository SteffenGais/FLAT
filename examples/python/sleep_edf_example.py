"""Download a public sleep recording and run sleep_example on it.

One night of the Sleep-EDF Database Expanded (PhysioNet, recording SC4001,
about 50 MB) is downloaded to examples/data, its hypnogram is converted into
a scoring text file, and the EEG channels Fpz-Cz and Pz-Oz are analyzed.

    python sleep_edf_example.py [DATA_FOLDER]

Data: Kemp et al. (2000), IEEE Trans Biomed Eng 47:1185-1194, and
Goldberger et al. (2000), Circulation 101:e215-e220 (license ODC-BY 1.0).
"""

import sys
from pathlib import Path
from urllib.request import urlretrieve

import matplotlib.pyplot as plt
import mne
import numpy as np

from sleep_example import sleep_example

DATA = Path(__file__).resolve().parents[1] / "data"
URL = "https://physionet-open.s3.amazonaws.com/sleep-edfx/1.0.0/sleep-cassette/"


def sleep_edf_example(data_folder=DATA):
    """Download recording SC4001 to data_folder, write its scoring file and analyze it."""
    data_folder = Path(data_folder)
    psg_file = data_folder / "SC4001E0-PSG.edf"
    hypnogram_file = data_folder / "SC4001EC-Hypnogram.edf"
    scoring_file = data_folder / "SC4001E0-scoring.txt"
    data_folder.mkdir(parents=True, exist_ok=True)
    for file in (psg_file, hypnogram_file):
        if not file.exists():
            print(f"Downloading {file.name}")
            part = file.with_name(file.name + ".part")
            urlretrieve(URL + file.name, part)
            part.rename(file)

    # One stage code per 30-s epoch of the recording: 0 W, 1-4 NREM, 5 REM,
    # 8 movement time, 9 not scored.
    psg = mne.io.read_raw_edf(psg_file, verbose="error")
    hypnogram_start = mne.io.read_raw_edf(hypnogram_file, verbose="error").info["meas_date"]
    offset = (hypnogram_start - psg.info["meas_date"]).total_seconds()
    annotations = mne.read_annotations(hypnogram_file)      # onsets from hypnogram_start
    codes = {"Sleep stage W": 0, "Sleep stage 1": 1, "Sleep stage 2": 2, "Sleep stage 3": 3,
             "Sleep stage 4": 4, "Sleep stage R": 5, "Movement time": 8}
    stage = np.full(psg.n_times // round(30 * psg.info["sfreq"]), 9)
    for onset, duration, text in zip(annotations.onset, annotations.duration,
                                     annotations.description):
        if text in codes:
            first = round((onset + offset) / 30)
            stage[first:first + round(duration / 30)] = codes[text]

    # The recording lasts almost a day. Analyze the night only, from 30 min
    # before sleep onset to 30 min after the final awakening.
    asleep = np.flatnonzero((stage >= 1) & (stage <= 5))
    stage[:max(asleep[0] - 60, 0)] = 9
    stage[asleep[-1] + 61:] = 9
    np.savetxt(scoring_file, stage, fmt="%d")

    # The EEG was high-pass filtered at 0.5 Hz: analyze from 1 Hz upwards.
    sleep_example(psg_file, scoring_file, channels=["EEG Fpz-Cz", "EEG Pz-Oz"],
                  frequency_range=(1, 45))


if __name__ == "__main__":
    sleep_edf_example(*sys.argv[1:2])
    plt.show()
