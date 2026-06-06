IAFout = readmatrix('IAFout_7ynew.xlsx');

% load IAFout.mat
cfg = [];


if ~isfield(cfg, 'age'),        cfg.age         = 84;    end; age         = cfg.age;
if ~isfield(cfg, 'bltype'),     cfg.bltype      = 'ecl'; end; bltype      = cfg.bltype;
if ~isfield(cfg, 'replacearg'), cfg.replacearg  = 0;     end; replacearg  = cfg.replacearg;
if ~isfield(cfg, 'plotarg'),    cfg.plotarg     = 0;     end; plotarg     = cfg.plotarg;
if ~isfield(cfg, 'segduration'),cfg.segduration = 2;     end; segduration = cfg.segduration;
if ~isfield(cfg, 'icapruning'), cfg.icapruning  = 1;     end; icapruning  = cfg.icapruning;

global EEG;

folderBase  = '/home/xielab/YiyiWang/FC analysis/Preprocessing';
folderout = fullfile(folderBase,'SegmentAverageFiles',filesep);
SegmentAverageFiles = fullfile(folderout,[num2str(age/12),'Y'],filesep);

for i = 8: length(IAFout)
    participantnumber = IAFout(i,1);

if participantnumber<10
    participantnumberstring = ['00' num2str(participantnumber)];
elseif participantnumber<100
    participantnumberstring = ['0' num2str(participantnumber)];
else
    participantnumberstring = num2str(participantnumber);
end

if age>12
    agestring  = [num2str(age./12) 'YF'];
    keystrings = {participantnumberstring,agestring,'.mat'};
    datapath   = fullfile(folderBase,'RawData',agestring,filesep);
else
    agestring = [num2str(age) 'mos'];
    keystrings = {participantnumberstring,'.mat'};
    datapath   = fullfile(folderBase,'RawData','infants',filesep);
end

Finaldataset = ['Expt 3 Participantnumber ' num2str(participantnumber) ' Age ' agestring ' Baseline_' bltype '_' num2str(segduration) 's_ica' num2str(icapruning) '.set'];


%load data
EEG = pop_loadset('filename',Finaldataset,'filepath',SegmentAverageFiles);
EEG = eeg_checkset(EEG);

folderBase  = '/home/xielab/YiyiWang/FC analysis/Preprocessing';
SegmentAverageFiles1 = fullfile(folderBase,'SegmentAverageFiles',filesep,'theta',filesep);
SegmentAverageFiles2 = fullfile(folderBase,'SegmentAverageFiles',filesep,'alpha',filesep);
SegmentAverageFiles3 = fullfile(folderBase,'SegmentAverageFiles',filesep,'beta',filesep);
SegmentAverageFiles4 = fullfile(folderBase,'SegmentAverageFiles',filesep,'gamma',filesep);

% filename
Finaldataset1 = ['theta ' num2str(participantnumber) ' Age ' agestring ' Baseline_' bltype '_' num2str(segduration) 's_ica' num2str(icapruning) '.set'];
Finaldataset2 = ['alpha ' num2str(participantnumber) ' Age ' agestring ' Baseline_' bltype '_' num2str(segduration) 's_ica' num2str(icapruning) '.set'];
Finaldataset3 = ['beta ' num2str(participantnumber) ' Age ' agestring ' Baseline_' bltype '_' num2str(segduration) 's_ica' num2str(icapruning) '.set'];
Finaldataset4 = ['gamma ' num2str(participantnumber) ' Age ' agestring ' Baseline_' bltype '_' num2str(segduration) 's_ica' num2str(icapruning) '.set'];



%% Bandpassing from here

% IAF
IAF = IAFout(find(IAFout(:,1)==participantnumber),2);
theta = [.4*IAF,.8*IAF];
alpha = [.8*IAF, 1.2*IAF];
beta = [1.2*IAF,30];
gamma = [30,50];

% 1. Filtering

% theta
EEG1  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff', theta, 'Design', 'butter', 'Filter', 'bandpass', 'Order',  8 ); 
EEG1  = pop_saveset(EEG1, 'filename',Finaldataset1,'filepath', SegmentAverageFiles1,'savemode','onefile');

% alpha
EEG2  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff', alpha, 'Design', 'butter', 'Filter', 'bandpass', 'Order',  8 ); 
EEG2  = pop_saveset(EEG2, 'filename',Finaldataset2,'filepath', SegmentAverageFiles2,'savemode','onefile');

% beta
EEG3  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff', beta, 'Design', 'butter', 'Filter', 'bandpass', 'Order',  8 ); 
EEG3  = pop_saveset(EEG3, 'filename',Finaldataset3,'filepath', SegmentAverageFiles3,'savemode','onefile');

% gamma
EEG4  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff', gamma, 'Design', 'butter', 'Filter', 'bandpass', 'Order',  8 ); 
EEG4  = pop_saveset(EEG4, 'filename',Finaldataset4,'filepath', SegmentAverageFiles4,'savemode','onefile');

end
