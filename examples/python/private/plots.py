"""Plots for the FLAT examples."""

import matplotlib.pyplot as plt
import numpy as np

from .recording import read_positions


def plot_spectra(freqs, power, baseline, flattened, channels, titles):
    """Line plots for one or two channels: power with its aperiodic baseline
    (top) and the flattened spectrum, power / baseline (bottom).

    power, baseline, flattened: conditions x channels x frequencies; one column
    of plots per condition.
    """
    n = len(titles)
    fig, axes = plt.subplots(2, n, figsize=(1.5 + 2.7 * n, 6.2), sharex=True,
                             sharey="row", squeeze=False, layout="constrained")
    for c, (top, bottom) in enumerate(axes.T):
        for k, name in enumerate(channels):
            top.plot(freqs, power[c, k], color=f"C{k}", linewidth=1.2, label=name)
            top.plot(freqs, baseline[c, k], "--", color=f"C{k}", linewidth=1)
            bottom.plot(freqs, flattened[c, k], color=f"C{k}", linewidth=1.2)
        bottom.axhline(1, color="0.5", linestyle=":", linewidth=1)
        top.set_yscale("log")
        top.set_title(titles[c])
        bottom.set_xlabel("Frequency (Hz)")
    ticks = [t for t in (0.5, 1, 2, 5, 10, 20, 40) if freqs[0] <= t <= freqs[-1]]
    for ax in axes.flat:
        ax.set_xscale("log")
        ax.set_xlim(freqs[0], freqs[-1])
        ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
        ax.set_xticks([], minor=True)
    axes[0, 0].set_ylabel("Power (µV²/Hz)")
    axes[1, 0].set_ylabel("Power / baseline")
    axes[0, 0].plot([], [], "k--", linewidth=1, label="Baseline")
    axes[0, 0].legend(loc="lower left")
    return fig


def plot_waterfall(freqs, power, flattened, channels, title):
    """Waterfall plots of power (left) and the flattened spectrum, power /
    baseline (right), for many channels. Channels with standard 10-05 names
    run from anterior (front) to posterior (back).

    power, flattened: channels x frequencies.
    """
    order = np.arange(len(channels))
    positions = read_positions()
    if all(name in positions for name in channels):
        order = np.argsort([-positions[name][1] for name in channels], kind="stable")
    keep = freqs <= 30
    fig = plt.figure(figsize=(16, 6.2))
    panels = ((10 * np.log10(power), "Power (dB)", "Power spectrum"),
              (flattened, "Power / baseline", "Flattened spectrum"))
    for k, (z, z_label, name) in enumerate(panels):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d", proj_type="ortho",
                             computed_zorder=False)
        draw_waterfall(ax, freqs[keep], z[order][:, keep], [channels[i] for i in order], z_label)
        ax.set_title(name, fontweight="bold")
    fig.suptitle(title)
    fig.subplots_adjust(left=0, right=0.96, bottom=0.03, top=0.93, wspace=0.04)
    return fig


def draw_waterfall(ax, freqs, z, channels, z_label):
    """Colored surface with a black line per channel, like MATLAB's surf with
    MeshStyle 'row'. It is drawn from back to front, so that front channels
    hide the channels behind them."""
    n = len(channels)
    rows = np.arange(1, n + 1)
    freq_grid, row_grid = np.meshgrid(freqs, rows)
    norm = plt.Normalize(np.nanmin(z), np.nanmax(z))
    layer = 0
    for k in range(n - 2, -1, -1):
        # The surface between channels k and k + 1, then the line of channel k + 1.
        ax.plot_surface(freq_grid[k:k + 2], row_grid[k:k + 2], z[k:k + 2], cmap="turbo",
                        norm=norm, rstride=1, cstride=1, linewidth=0,
                        antialiased=False, shade=False, zorder=layer)
        ax.plot(freqs, np.full(freqs.size, rows[k + 1]), z[k + 1], color="k",
                linewidth=0.4, zorder=layer + 1)
        layer += 2
    ax.plot(freqs, np.full(freqs.size, rows[0]), z[0], color="k", linewidth=0.4, zorder=layer)
    ax.view_init(elev=50, azim=-81)
    ax.set_box_aspect((2, 1.3, 0.9), zoom=1.15)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((1, 1, 1, 1))
    step = int(np.ceil(n / 12))
    ax.set_xlim(0, 30)
    ax.set_ylim(1, n)
    ax.set_yticks(rows[::step], channels[::step])
    ax.set_xlabel("Frequency (Hz)")
    ax.set_zlabel(z_label)
