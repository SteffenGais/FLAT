function [power, freqs] = welch_psd(x, fs, seconds)
%WELCH_PSD Power spectral density by Welch's method.
%   x: channels x samples. Hann windows of the given length in seconds,
%   50% overlap, linear detrending of each window. Returns the one-sided
%   density (channels x frequencies), e.g. in µV^2/Hz for data in µV.
%   Same as scipy.signal.welch(x, fs, nperseg=n, detrend='linear').
n = round(seconds * fs);
starts = 1 : n - floor(n / 2) : size(x, 2) - n + 1;
if isempty(starts)
    error('The data are shorter than one %g-s window.', seconds);
end
taper = 0.5 - 0.5 * cos(2 * pi * (0:n-1)' / n);     % periodic Hann window
power = 0;
for s = starts
    segment = detrend(x(:, s:s+n-1)') .* taper;     % samples x channels
    power = power + abs(fft(segment)).^2;
end
power = power(1:floor(n/2) + 1, :)' / (numel(starts) * fs * sum(taper.^2));
last = size(power, 2) - (mod(n, 2) == 0);           % Nyquist bin counts once
power(:, 2:last) = 2 * power(:, 2:last);
freqs = (0:size(power, 2) - 1) * fs / n;
end
