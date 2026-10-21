"""Power spectra for the FLAT examples."""

import numpy as np


def log_cells(freqs, power):
    """Average power within 1/32-octave cells.

    Groups the frequency bins by their 1/32-octave cell and averages frequency
    and power (linear units) within each cell. Where a cell holds a single bin
    (low frequencies), the bin is kept. The result has constant relative
    resolution: full resolution at low frequencies, many averaged bins at high
    frequencies. power: (..., frequencies).
    """
    positive = freqs > 0
    f, p = freqs[positive], power[..., positive]
    cells = np.round(32 * np.log2(f)).astype(int)
    _, first, counts = np.unique(cells, return_index=True, return_counts=True)
    return np.add.reduceat(f, first) / counts, np.add.reduceat(p, first, axis=-1) / counts
