# FLAT method
[Overview](../README.md) · [MATLAB](../README_matlab.md) · [Python](../README_python.md)

FLAT (Flattening spectra using the Local Aperiodic Trend) separates a power
spectrum into a smooth aperiodic background and the periodic activity on top
of it.

## Decomposition

For a power spectrum $P(f)$, FLAT estimates the aperiodic background $B(f)$,
the *baseline*, and returns the flattened spectrum $R(f) = P(f) / B(f)$.
The baseline has the units of $P$; $R$ is dimensionless and equals 1 where
there is only background. In log power, $\log P = \log B + \log R$.
The flattened spectrum shows every deviation from the background: periodic
activity as peaks, together with the noise of the estimate. No threshold or
peak detection decides what is shown.

The input must be linear power or power spectral density, not amplitude or
decibels. FLAT analyzes the supplied frequencies within **0.3–45 Hz**; all
other frequencies are NaN in the outputs. Multiplying the input by a constant
(e.g. changing units) scales the baseline and leaves $R$ and the slope unchanged.

## Algorithm

1. Keep the frequencies within 0.3–45 Hz. Bins inside declared notch bands
   are bridged by linear interpolation of log power over log frequency from
   their unaffected neighbors. They support the estimate but are NaN in both
   outputs.
2. With $x = \log_2 f$ and $y = \log_{10} P$, interpolate $y$ onto a grid with
   1/32-octave spacing, starting at the first frequency.
3. Fit a Theil–Sen line to the outer 32 grid points (one octave) at each end:
   its slope is the median of all pairwise slopes, its intercept the median
   residual. Extend the spectrum along these lines by 65 points on both sides.
4. Apply a centered moving median over 65 points (two octaves), followed by a
   moving mean over 65 points. The gradient of this smooth curve is the local
   slope $s_i$.
5. In each window, take the median of $y_j - s_i \cdot (j - i)$, i.e. of the
   values after removing the local slope. Smooth these medians with a moving
   mean over 65 points. This is the log baseline.
6. Remove the extensions, interpolate the log baseline back to the original
   frequencies and exponentiate it to obtain $B$. Return $P / B$.
7. The slope summary regresses $\log_{10} B$ on $\log_{10} f$ at 116 log-spaced
   frequencies from **2 to 24 Hz**. A $1/f^2$ background gives −2. For a power
   law $1/f^\chi$, the slope is $-\chi$; the aperiodic exponent reported by
   specparam (FOOOF) in its fixed mode therefore corresponds to −slope.

## Choice of the fixed constants

- **Two-octave window.** FLAT counts activity as periodic if its spectral peak
  is narrower than one octave at half maximum. A peak one octave wide
  corresponds to a quality factor Q ≈ 1.4: such a rhythm loses its phase
  within about half a cycle, so wider activity is not periodic in any useful
  sense. A median stays within the background values as long as a peak
  covers less than half of its window, but a strong peak rises above the
  noise over two to three times its half-maximum width. The two-octave window
  therefore keeps narrow rhythms such as sleep spindles and alpha (about
  0.2–0.4 octave at half maximum) largely out of the baseline, while wider
  peaks are partly absorbed (see Limits). A narrower window would follow the
  background more closely but would be lifted by more peaks. Consequently,
  broad activity, such as most delta activity in deep sleep, is part of the
  baseline, whereas narrow rhythms such as alpha, sleep spindles or a
  slow-oscillation peak near 1 Hz remain in the flattened spectrum.
- **1/32-octave grid.** Each window holds 65 points, enough for a stable median.
  The grid is finer than the FFT bins at low frequencies and comparable to
  them around 10 Hz.
- **Edge fit over one octave.** Half a window: long enough for a robust line,
  short enough to describe only the edge of the spectrum.
- **Analysis range 0.3–45 Hz.** The range most EEG and MEG recordings measure
  reliably: above slow drifts and typical high-pass filters, below line noise
  at 50 or 60 Hz.
- **Slope interval 2–24 Hz.** Above the high-pass and slow-oscillation region
  (the baseline at 2 Hz uses data down to 1 Hz) and below most muscle activity
  and the 35-Hz low-pass filter common in sleep recordings. The two-octave
  window around 24 Hz ends at 48 Hz, below 50-Hz line noise.

## Properties

- **No model and nothing to tune.** FLAT fits no peaks and assumes no shape
  for the background as a whole. There are no settings for peak number, peak
  width, thresholds or the aperiodic model; the same call works for every
  recording.
- **Local background.** Within a two-octave window, the background is treated
  as a straight line in log-log coordinates. The baseline therefore follows
  backgrounds that bend, such as spectra with a knee, rather than forcing a
  single power law onto the whole spectrum.
- **Robust to narrow peaks.** A median depends on how much of the window a
  peak covers, not on how high it is, so narrow peaks hardly move the
  baseline. Removing the local slope before taking the median (step 5) keeps
  a peak on a sloping background from lifting the median to a higher part of
  the slope.
- **Full window at the edges.** The robust continuation (step 3) keeps the
  full two-octave window at both ends of the range instead of letting it
  shrink.
- **Equal weight per octave.** On the log-frequency grid, high frequencies do
  not dominate the estimate merely because a linear FFT grid has more bins
  there. The same holds for the slope.
- **Any power spectrum.** FLAT needs only the spectrum (e.g. from Welch's
  method, multitapers or wavelets) within the analyzed range, not the time
  series or data beyond the range.

## Input spectra

FLAT works on a log-frequency axis, so the input spectrum should have
constant relative resolution. Short windows give few frequency bins at low
frequencies (with 5-s windows, the octave from 0.3 to 0.6 Hz contains only
two bins); long windows give noisy high frequencies. Both are avoided by using
windows as long as the recording allows while still averaging enough of them
(about 30), e.g. one 30-s Hann window per sleep epoch, or 4-s windows with
50% overlap for one minute of data. The periodograms are averaged in linear
power, and power is then averaged within 1/32-octave cells: at low
frequencies every bin is kept, at high frequencies many bins are averaged.
In sleep recordings, this gave the same results as 5-s windows above 3 Hz,
resolved the range below 1 Hz much better and was faster to compute. The
READMEs show the code.

## Reliable frequency range

Filters bend the spectrum near their cutoff frequencies, and the two-octave
window carries this bend about one octave into the estimate. Before calling
FLAT, keep only the frequencies the recording measures reliably, with a margin
to the filter edges. Trim power and frequencies together:

```matlab
keep = freqs >= 1 & freqs <= 40;
[flattened, baseline, slope] = flatten_spectrum(power(:, keep), freqs(keep));
```

```python
keep = (freqs >= 1) & (freqs <= 40)
flattened, baseline, slope = flatten_spectrum(power[..., keep], freqs[keep])
```

The slope needs the whole 2–24 Hz interval. It is NaN if the frequencies
start above 2 Hz or end below 24 Hz, or if a notch band overlaps this
interval. The baseline and the flattened spectrum are still estimated.

## Limits

- **Too little averaging.** The median of a noisy spectrum lies below its
  mean, which lowers the baseline and raises the flattened spectrum. For an
  average of $K$ independent periodograms, the baseline is about $1/(3K)$ too
  low: about 3% for $K = 10$ and 0.3% for $K = 100$. Average many windows in
  linear power (e.g. Welch's method) before flattening.
- **Broad or neighboring peaks.** Peaks are partly absorbed into the baseline
  when they, or several neighboring peaks together, cover a large part of the
  two-octave window. In simulations, a peak five times the background came
  out about 7% too small at a half-maximum width of 0.25 octave and about 23%
  too small at 0.5 octave. Compare power and baseline when interpreting broad
  elevations.
- **Peaks at the ends of the range.** The edge continuation is fitted to the
  outermost octave, so a strong peak there tilts it and is partly absorbed.
  In simulations, a slow-oscillation peak at 0.9 Hz was unaffected when the
  analysis started at 0.6 Hz or lower, but came out 30–65% too small when it
  started at 0.7–0.8 Hz. Start the analysis about one octave below the lowest
  peak of interest, e.g. at 0.3–0.5 Hz for slow oscillations near 1 Hz. This
  needs a correspondingly low high-pass filter and long windows. 
- **Strong curvature.** Curvature within two octaves is approximated by
  local straight lines, and a single slope cannot describe a curved
  background. Plot the baseline when curvature matters.
