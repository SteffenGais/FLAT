function sleep_edf_example(data_folder)
%SLEEP_EDF_EXAMPLE Download a public sleep recording and run SLEEP_EXAMPLE.
%   SLEEP_EDF_EXAMPLE downloads one night of the Sleep-EDF Database Expanded
%   (PhysioNet, recording SC4001, about 50 MB) to examples/data, converts its
%   hypnogram into a scoring text file and runs SLEEP_EXAMPLE on the EEG
%   channels Fpz-Cz and Pz-Oz.
%
%   SLEEP_EDF_EXAMPLE(DATA_FOLDER) downloads to another folder.
%
%   Data: Kemp et al. (2000), IEEE Trans Biomed Eng 47:1185-1194, and
%   Goldberger et al. (2000), Circulation 101:e215-e220 (license ODC-BY 1.0).
arguments
    data_folder {mustBeTextScalar} = fullfile(fileparts(mfilename('fullpath')), '..', 'data')
end
url = 'https://physionet-open.s3.amazonaws.com/sleep-edfx/1.0.0/sleep-cassette/';
psg_file = fullfile(data_folder, 'SC4001E0-PSG.edf');
hypnogram_file = fullfile(data_folder, 'SC4001EC-Hypnogram.edf');
scoring_file = fullfile(data_folder, 'SC4001E0-scoring.txt');
if ~isfolder(data_folder)
    mkdir(data_folder);
end
for file = string({psg_file, hypnogram_file})
    if ~isfile(file)
        [~, name, extension] = fileparts(file);
        fprintf('Downloading %s%s\n', name, extension);
        websave(file + ".part", url + name + extension);
        movefile(file + ".part", file);
    end
end

% One stage code per 30-s epoch of the recording: 0 W, 1-4 NREM, 5 REM,
% 8 movement time, 9 not scored.
psg = read_eeg(psg_file);
hypnogram = read_eeg(hypnogram_file);
labels = ["Sleep stage W", "Sleep stage 1", "Sleep stage 2", "Sleep stage 3", ...
    "Sleep stage 4", "Sleep stage R", "Movement time"];
codes = [0 1 2 3 4 5 8];
n_epochs = floor(psg.n_samples / psg.fs / 30);
stage = 9 * ones(n_epochs, 1);
offset = seconds(hypnogram.start - psg.start);
for k = 1:height(hypnogram.annotations)
    a = hypnogram.annotations(k, :);
    [known, c] = ismember(a.text, labels);
    if known
        first = round((a.onset + offset) / 30) + 1;
        last = min(first + round(a.duration / 30) - 1, n_epochs);
        stage(first:last) = codes(c);
    end
end

% The recording lasts almost a day. Analyze the night only, from 30 min
% before sleep onset to 30 min after the final awakening.
asleep = find(stage >= 1 & stage <= 5);
stage([1:asleep(1) - 61, asleep(end) + 61:end]) = 9;
writematrix(stage, scoring_file, 'FileType', 'text');

% The EEG was high-pass filtered at 0.5 Hz: analyze from 1 Hz upwards.
sleep_example(psg_file, scoring_file, ...
    'Channels', ["EEG Fpz-Cz", "EEG Pz-Oz"], 'FrequencyRange', [1 45]);
end
