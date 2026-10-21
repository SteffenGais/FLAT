"""FLAT: Flattening spectra using the Local Aperiodic Trend."""

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view


def flatten_spectrum(power, freqs, notch_hz=None):
    """Return flattened power / baseline, aperiodic baseline, and 2–24 Hz slope.

    Parameters
    ----------
    power : array_like, shape (..., n_frequencies)
        Linear power or power spectral density; frequencies are the last axis.
        An entirely NaN spectrum is allowed and returns entirely NaN outputs.
    freqs : array_like, shape (n_frequencies,)
        Finite, nonnegative, strictly increasing frequencies in Hz.
    notch_hz : None, pair, or array_like of shape (n_bands, 2), optional
        Excluded intervals in Hz, e.g. (32, 33). None (default) or [] means
        no notch. Give the actual affected interval, not just its center.

    Returns
    -------
    flattened, baseline : ndarray, same shape as power
        Dimensionless power / baseline and baseline in the input power units.
        Values outside 0.3–45 Hz and inside notch bands are NaN.
    slope : ndarray, shape power.shape[:-1]
        Signed least-squares slope of log10(baseline) versus log10(frequency),
        with equal weight per log-frequency interval over 2–24 Hz. A 1/f**2
        baseline has slope -2. NaN if this full interval is unavailable or
        intersects a notch. For one input spectrum, this is a scalar array.

    Notes
    -----
    The two-octave window and 1/32-octave grid are fixed method constants.
    Supply averaged spectra from a reliable passband; FLAT does not detect
    filters or reject artifacts. See docs/method.md for assumptions and limits.
    """
    power = np.asarray(power, dtype=float)
    freqs = np.asarray(freqs, dtype=float)
    if (freqs.ndim != 1 or freqs.size < 2 or not np.isfinite(freqs).all()
            or np.any(freqs < 0) or np.any(np.diff(freqs) <= 0)):
        raise ValueError("freqs must be finite, nonnegative, and strictly increasing.")
    if power.ndim < 1 or power.shape[-1] != freqs.size or power.size == 0:
        raise ValueError("The last power axis must match freqs; power must be nonempty.")

    bands = np.asarray([] if notch_hz is None else notch_hz, dtype=float)
    if bands.size:
        bands = np.atleast_2d(bands)
        if (bands.ndim != 2 or bands.shape[1] != 2 or not np.isfinite(bands).all()
                or np.any(bands[:, 0] <= 0) or np.any(bands[:, 1] <= bands[:, 0])
                or np.any(bands[1:, 0] <= bands[:-1, 1])):
            raise ValueError("Notch bands must be positive, ordered, disjoint [low, high] pairs.")
    else:
        bands = np.empty((0, 2))

    # Work only in the fixed EEG/MEG analysis range. DC is never logarithmed.
    selected = (freqs >= 0.3) & (freqs <= 45)
    f = freqs[selected]
    if f.size < 2:
        raise ValueError("At least two frequency bins are required within 0.3–45 Hz.")
    x = np.log2(f)
    grid = np.arange(x[0], x[-1] + 1e-9, 1 / 32)
    if grid.size < 2:
        raise ValueError("The analyzed frequency span must cover at least 1/32 octave.")
    window = 65                         # Two octaves: 64 intervals, 65 centered points.
    half = window // 2
    n_fit = min(half, grid.size)          # Fit the outer one octave at each edge.
    offsets = np.arange(-half, half + 1)
    i, j = np.triu_indices(n_fit, 1)
    edge_x = np.arange(n_fit)
    masked = np.zeros(f.size, dtype=bool)
    for low, high in bands:
        masked |= (f >= low) & (f <= high)
    if masked.any() and (masked[0] or masked[-1]):
        raise ValueError("Each analyzed notch needs a frequency bin on both sides.")

    rows = power.reshape(-1, freqs.size)
    flattened = np.full_like(rows, np.nan)
    baseline = np.full_like(rows, np.nan)
    slopes = np.full(rows.shape[0], np.nan)
    for r, spectrum in enumerate(rows):
        if np.isnan(spectrum).all():
            continue
        p = spectrum[selected]
        if not np.isfinite(p[~masked]).all() or np.any(p[~masked] <= 0):
            raise ValueError("Analyzed power outside notches must be finite and positive.")

        # Bridge excluded bins from unaffected neighbors in log-log coordinates.
        # They support the baseline calculation but are never reported as data.
        logp = np.interp(x, x[~masked], np.log10(p[~masked]))
        y = np.interp(grid, x, logp)

        # Robust linear continuation prevents shrinking-window bias at the edges.
        # Theil–Sen slope = median of all pairwise slopes; intercept = median residual.
        extensions = []
        for edge, query in ((y[:n_fit], np.arange(-window, 0)),
                            (y[-n_fit:], np.arange(n_fit, n_fit + window))):
            b = np.median((edge[j] - edge[i]) / (j - i))
            a = np.median(edge - b * edge_x)
            extensions.append(a + b * query)
        y = np.concatenate((extensions[0], y, extensions[1]))

        # First median-then-mean pass estimates a smooth local slope. NaN padding
        # reproduces centered, shrinking windows at the far extrapolated ends.
        neighborhoods = sliding_window_view(np.pad(y, half, constant_values=np.nan), window)
        first = np.nanmedian(neighborhoods, axis=-1)
        first = np.nanmean(sliding_window_view(
            np.pad(first, half, constant_values=np.nan), window), axis=-1)
        local_slope = np.gradient(first)

        # Remove the local slope BEFORE taking the median: on a sloped spectrum,
        # a positive peak can otherwise shift the median into a higher background.
        detrended = neighborhoods - local_slope[:, None] * offsets
        trend = np.nanmedian(detrended, axis=-1)
        trend = np.nanmean(sliding_window_view(
            np.pad(trend, half, constant_values=np.nan), window), axis=-1)
        trend = trend[window:window + grid.size]

        # Interpolate back. The final original bin may lie just beyond the grid;
        # continue its last segment instead of silently clamping the endpoint.
        logbase = np.interp(x, grid, trend)
        beyond = x > grid[-1]
        logbase[beyond] = trend[-1] + (x[beyond] - grid[-1]) * (
            trend[-1] - trend[-2]) / (grid[-1] - grid[-2])
        b = 10 ** logbase
        b[masked] = np.nan
        baseline[r, selected] = b
        flattened[r, selected] = p / b

        # Equal log-frequency weighting makes the summary independent of the
        # density of linear FFT bins. Do not shorten an unavailable 2–24 Hz fit.
        if f[0] <= 2 and f[-1] >= 24 and not any(
                low <= 24 and high >= 2 for low, high in bands):
            fit_x = np.linspace(1, np.log2(24), 116)
            fit_y = np.interp(fit_x, x, logbase)
            centered = fit_x - fit_x.mean()
            slopes[r] = np.dot(centered, fit_y) / np.dot(centered, centered) / np.log10(2)

    return (flattened.reshape(power.shape), baseline.reshape(power.shape),
            slopes.reshape(power.shape[:-1]))
