function [freqs, power] = log_cells(freqs, power)
%LOG_CELLS Average power within 1/32-octave cells.
%   [freqs, power] = log_cells(freqs, power) groups the frequency bins by
%   their 1/32-octave cell and averages frequency and power (linear units)
%   within each cell. Where a cell holds a single bin (low frequencies), the
%   bin is kept. The result has constant relative resolution: full resolution
%   at low frequencies, many averaged bins at high frequencies.
%   power: channels x frequencies.
positive = freqs > 0;
freqs = freqs(positive);
power = power(:, positive);
[~, ~, cells] = unique(round(32 * log2(freqs)));
members = sparse(1:numel(cells), cells, 1);    % bins x cells
counts = full(sum(members, 1));
freqs = (freqs * members) ./ counts;
power = full(power * members) ./ counts;
end
