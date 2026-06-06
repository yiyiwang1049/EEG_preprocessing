function BPCD = SourceAnalysis2(cfg,participantnumber)


if ~isfield(cfg, 'age'),        cfg.age         = 84;    end; age         = cfg.age;
if ~isfield(cfg, 'bltype'),     cfg.bltype      = 'ecl'; end; bltype      = cfg.bltype;
if ~isfield(cfg, 'replacearg'), cfg.replacearg  = 0;     end; replacearg  = cfg.replacearg;
if ~isfield(cfg, 'plotarg'),    cfg.plotarg     = 0;     end; plotarg     = cfg.plotarg;
if ~isfield(cfg, 'segduration'),cfg.segduration = 2;     end; segduration = cfg.segduration;
if ~isfield(cfg, 'icapruning'), cfg.icapruning  = 1;     end; icapruning  = cfg.icapruning;
if ~isfield(cfg, 'band'),       cfg.band       = 'theta';end; band        = cfg.band;

if age>12
    agestring  = [num2str(age./12) 'YF'];
else
    agestring = [num2str(age) 'mos'];
end

Finaldataset = [band ' ' num2str(participantnumber) ' Age ' agestring ' Baseline_' bltype '_' num2str(segduration) 's_ica' num2str(icapruning) '.set'];
folderBase  = '/home/xielab/YiyiWang/FC analysis/Preprocessing';
folderout = fullfile(folderBase,'SegmentAverageFiles',filesep);
SegmentAverageFiles = fullfile(folderout,band,filesep);


%load data
EEG = pop_loadset('filename',Finaldataset,'filepath',SegmentAverageFiles);
EEG = eeg_checkset(EEG);

% combine all trials
a = EEG.data(:,:,1);
if size(EEG.data,3)>60 %max 3YF trials = 60 
    for i = 60:size(EEG.data,3)-1
        a = [a,EEG.data(:,:,i+1)];
    end
else
    BPCD = 0;
    return;
end

EEGcombine = a;

%  Hilbert transformed EEG data
EEG_hilbert = hilbert(EEGcombine);

% ource-space analytical time series
load('kernel_7Y.mat');
%kernel_5Y=kernel_5Y.ImagingKernel;
EEG_source = kernel_7Y*EEG_hilbert;

BPCD = EEG_source;
