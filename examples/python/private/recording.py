"""Reading recordings for the FLAT examples."""

from pathlib import Path

import mne
import numpy as np

POSITIONS = Path(__file__).with_name("standard_1005.tsv")


def read_positions():
    """Standard 10-05 electrode positions in mm (x right, y anterior, z up)."""
    lines = POSITIONS.read_text(encoding="utf-8").splitlines()[1:]
    return {name: np.array([float(x), float(y), float(z)])
            for name, x, y, z in (line.split("\t") for line in lines)}


def read_eeg(eeg_file, channels=None):
    """Open a recording with MNE, e.g. BrainVision (.vhdr) or EDF, without loading it.

    Channel names get the standard 10-05 spelling where possible, e.g. 'Fc5.'
    becomes 'FC5'. channels: names (not case-sensitive); None selects all.
    """
    raw = mne.io.read_raw(eeg_file, verbose="error")
    standard = {name.lower(): name for name in read_positions()}
    raw.rename_channels(lambda name: standard.get(name.replace(".", "").lower(), name))
    if channels is not None:
        names = {name.lower(): name for name in raw.ch_names}
        missing = [channel for channel in channels if channel.lower() not in names]
        if missing:
            raise ValueError(f"Channel(s) not found: {', '.join(missing)}")
        raw.pick([names[channel.lower()] for channel in channels])
    return raw
