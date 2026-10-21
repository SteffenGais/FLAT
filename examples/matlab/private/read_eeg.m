function rec = read_eeg(eeg_file, channels)
%READ_EEG Open a BrainVision (.vhdr) or EDF (.edf) recording.
%   rec = read_eeg(eeg_file, channels) returns a struct with
%     labels       channel names (standard 10-05 spelling where possible)
%     fs           sampling rate in Hz
%     n_samples    number of samples per channel
%     read(a, b)   samples a to b of the selected channels in µV
%                  (channels x samples)
%   EDF files also give
%     start        start of the recording (datetime)
%     annotations  EDF+ annotations: onset and duration in s, text
%   channels: "all" (default), names, or channel numbers.
%   EDF signals with a lower sampling rate (e.g. respiration) are skipped.
if nargin < 2
    channels = "all";
end
[~, ~, extension] = fileparts(eeg_file);
switch lower(extension)
    case '.vhdr'
        rec = read_brainvision(eeg_file);
    case '.edf'
        rec = read_edf(eeg_file);
    otherwise
        error('Please use a BrainVision (.vhdr) or EDF (.edf) file.');
end

% Standard 10-05 spelling, e.g. 'Fc5.' becomes 'FC5'.
positions = readtable(fullfile(fileparts(mfilename('fullpath')), 'standard_1005.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
[known, k] = ismember(lower(erase(rec.labels, '.')), lower(positions.name));
rec.labels(known) = positions.name(k(known));

if isnumeric(channels)
    picks = channels(:)';
elseif isequal(lower(string(channels)), "all")
    picks = 1:numel(rec.labels);
else
    requested = string(channels);
    [found, picks] = ismember(lower(requested), lower(string(rec.labels)));
    if ~all(found)
        error('Channel(s) not found: %s', strjoin(requested(~found), ', '));
    end
end
rec.labels = rec.labels(picks);
read_channels = rec.read;
rec.read = @(first, last) read_channels(first, last, picks);
end


function rec = read_brainvision(vhdr_file)
header = fileread(vhdr_file);
field = @(key) strtrim(char(regexp(header, ['(?m)^' key '=([^\r\n]*)'], 'tokens', 'once')));
layout.n_channels = str2double(field('NumberOfChannels'));
rec.fs = 1e6 / str2double(field('SamplingInterval'));
switch field('BinaryFormat')
    case 'INT_16'
        layout.precision = 'int16=>double';
        layout.bytes = 2;
    case 'IEEE_FLOAT_32'
        layout.precision = 'single=>double';
        layout.bytes = 4;
    otherwise
        error('Unsupported BrainVision binary format: %s', field('BinaryFormat'));
end
layout.file = fullfile(fileparts(vhdr_file), field('DataFile'));
info = dir(layout.file);
rec.n_samples = floor(info.bytes / layout.bytes / layout.n_channels);
layout.n_samples = rec.n_samples;
layout.multiplexed = ~strcmpi(field('DataOrientation'), 'VECTORIZED');

% Ch<number>=<name>,<reference>,<resolution>,<unit>
rec.labels = cell(1, layout.n_channels);
layout.scale = ones(layout.n_channels, 1);
for c = 1:layout.n_channels
    parts = split(string(field(sprintf('Ch%d', c))), ',');
    rec.labels{c} = char(strrep(parts(1), '\1', ','));
    if numel(parts) >= 3 && strlength(parts(3)) > 0
        layout.scale(c) = str2double(parts(3));
    end
    if numel(parts) >= 4
        layout.scale(c) = layout.scale(c) * to_microvolt(parts(4));
    end
end
rec.read = @(first, last, picks) read_brainvision_samples(layout, first, last, picks);
end


function x = read_brainvision_samples(layout, first, last, picks)
fid = fopen(layout.file, 'r', 'ieee-le');
closer = onCleanup(@() fclose(fid));
n = last - first + 1;
if layout.multiplexed
    fseek(fid, (first - 1) * layout.n_channels * layout.bytes, 'bof');
    x = fread(fid, [layout.n_channels, n], layout.precision);
    x = x(picks, :);
else
    x = zeros(numel(picks), n);
    for k = 1:numel(picks)
        fseek(fid, ((picks(k) - 1) * layout.n_samples + first - 1) * layout.bytes, 'bof');
        x(k, :) = fread(fid, [1, n], layout.precision);
    end
end
x = x .* layout.scale(picks);
end


function rec = read_edf(edf_file)
fid = fopen(edf_file, 'r', 'ieee-le');
closer = onCleanup(@() fclose(fid));
header = fread(fid, [1, 256], '*char');
n_signals = str2double(header(253:256));
record_seconds = str2double(header(245:252));
rec.start = datetime([header(169:176) ' ' header(177:184)], ...
    'InputFormat', 'dd.MM.yy HH.mm.ss', 'PivotYear', 1985);

% Signal header fields: label, transducer, unit, physical min, physical
% max, digital min, digital max, prefiltering, samples per record, reserved.
widths = [16 80 8 8 8 8 8 80 8 32];
signal_header = fread(fid, [1, 256 * n_signals], '*char');
fields = cell(1, numel(widths));
position = 0;
for k = 1:numel(widths)
    block = signal_header(position + (1:widths(k) * n_signals));
    fields{k} = strtrim(cellstr(reshape(block, widths(k), n_signals)'));
    position = position + widths(k) * n_signals;
end
labels = fields{1};
physical = str2double([fields{4}, fields{5}]);
digital = str2double([fields{6}, fields{7}]);
layout.per_record = str2double(fields{9});
layout.file = edf_file;
layout.header_bytes = 256 * (n_signals + 1);
layout.record_values = sum(layout.per_record);
layout.offset = [0; cumsum(layout.per_record(1:end-1))];

% Value in µV = gain * stored value + shift.
scale = cellfun(@to_microvolt, fields{3});
layout.gain = scale .* diff(physical, 1, 2) ./ diff(digital, 1, 2);
layout.shift = scale .* physical(:, 1) - layout.gain .* digital(:, 1);
info = dir(edf_file);
n_records = floor((info.bytes - layout.header_bytes) / (2 * layout.record_values));

% Data channels: all ordinary signals with the highest sampling rate.
is_annotation = strcmp(labels, 'EDF Annotations');
rate = layout.per_record .* ~is_annotation;
layout.data = find(rate > 0 & rate == max(rate));
rec.labels = labels(layout.data)';
if isempty(layout.data)
    rec.fs = NaN;
    rec.n_samples = 0;
else
    rec.fs = layout.per_record(layout.data(1)) / record_seconds;
    rec.n_samples = n_records * layout.per_record(layout.data(1));
end
rec.read = @(first, last, picks) read_edf_samples(layout, first, last, picks);

% EDF+ annotations: +onset[\x15duration]\x14text\x14[text\x14...]\x00
rec.annotations = table('Size', [0 3], 'VariableTypes', {'double', 'double', 'string'}, ...
    'VariableNames', {'onset', 'duration', 'text'});
for s = find(is_annotation)'
    fseek(fid, layout.header_bytes, 'bof');
    values = fread(fid, [layout.record_values, n_records], 'int16=>int16');
    bytes = typecast(reshape(values(layout.offset(s) + (1:layout.per_record(s)), :), 1, []), 'uint8');
    for tal = split(string(char(bytes)), char(0))'
        parts = split(tal, char(20));
        if numel(parts) < 3
            continue
        end
        time = [str2double(split(parts(1), char(21))); 0];
        for annotation = parts(2:end-1)'
            if strlength(annotation) > 0
                rec.annotations(end+1, :) = {time(1), time(2), annotation};
            end
        end
    end
end
end


function x = read_edf_samples(layout, first, last, picks)
signals = layout.data(picks);
n = layout.per_record(signals(1));
records = floor((first - 1) / n) : floor((last - 1) / n);
fid = fopen(layout.file, 'r', 'ieee-le');
closer = onCleanup(@() fclose(fid));
fseek(fid, layout.header_bytes + records(1) * layout.record_values * 2, 'bof');
values = fread(fid, [layout.record_values, numel(records)], 'int16=>double');
x = zeros(numel(signals), n * numel(records));
for k = 1:numel(signals)
    x(k, :) = reshape(values(layout.offset(signals(k)) + (1:n), :), 1, []);
end
x = x(:, first - records(1) * n : last - records(1) * n);
x = x .* layout.gain(signals) + layout.shift(signals);
end


function factor = to_microvolt(unit)
% Factor from a voltage unit to µV. Anything else, e.g. 'uV' or 'µV' in
% any text encoding, is taken as µV.
switch strtrim(char(unit))
    case 'V'
        factor = 1e6;
    case 'mV'
        factor = 1e3;
    case 'nV'
        factor = 1e-3;
    otherwise
        factor = 1;
end
end
