% This program will read in the .mat file for the EmoKidsFC project and perform data preprocessing
%----------------------------------------------------------------------------------------
% % For Baseline analysis:
% 1) read the raw Matlab .mat file
% 2) extract the baseline section 
% 3) preprocessing and save the processed data
%-------------------------------------------------------------------------------------------------
% Note for session 3: Events are "words" and marked at the beginning of the seg (baseline) not the stimulus;
function [] = ReadinDataAndPreprocessing_baseline(cfg,participantnumber)
if ~isfield(cfg, 'age'),        cfg.age         = 84;    end; age         = cfg.age;
if ~isfield(cfg, 'bltype'),     cfg.bltype      = 'eop'; end; bltype      = cfg.bltype;
if ~isfield(cfg, 'replacearg'), cfg.replacearg  = 0;     end; replacearg  = cfg.replacearg;
if ~isfield(cfg, 'plotarg'),    cfg.plotarg     = 0;     end; plotarg     = cfg.plotarg;
if ~isfield(cfg, 'segduration'),cfg.segduration = 2;     end; segduration = cfg.segduration;
if ~isfield(cfg, 'icapruning'), cfg.icapruning  = 1;     end; icapruning  = cfg.icapruning;

%eeglab; close; ft_defaults;
global EEG;

folderBase  = '/home/xielab/YiyiWang/FC analysis/Preprocessing';
if cfg.age == 36
SegmentAverageFiles0 = fullfile(folderBase,'SegmentAverageFiles',filesep);
SegmentAverageFiles = fullfile(SegmentAverageFiles0,'3Y',filesep);
end
if cfg.age == 60
SegmentAverageFiles0 = fullfile(folderBase,'SegmentAverageFiles',filesep);
SegmentAverageFiles = fullfile(SegmentAverageFiles0,'5Y',filesep);
end
if cfg.age == 84
SegmentAverageFiles0 = fullfile(folderBase,'SegmentAverageFiles',filesep);
SegmentAverageFiles = fullfile(SegmentAverageFiles0,'7Y',filesep);
end

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

%% check if the output EEGLAB data already exists
Finaldataset = ['Expt 3 Participantnumber ' num2str(participantnumber) ' Age ' agestring ' Baseline_' bltype '_' num2str(segduration) 's_ica' num2str(icapruning) '.set'];
if exist([SegmentAverageFiles Finaldataset]) && replacearg == 0
    disp([Finaldataset ' already exists. Change cfg.replacearg to 1 for replacement.']);
    return
end

%% load the data
filename   = find_filename(keystrings,datapath);
% read in the egi .mat file
fieldname  = strrep(filename,'.mat','');
srate      = 500.00;
% NOTE: WX changed the default in ‘readegilocs.m’
% pop_importegimat calls readegilocs.m which loads the GSN129.sfg as default;
% WX changed the default in readegilocs.m for EGI data to GSN-HydroCel-129.sfp
EEG = pop_importegimat([datapath filename], srate, 0, fieldname);
%EEG = pop_importegimat([datapath filename], srate, 0);
% disp(EEG.chaninfo.filename);

% load channel
    channel_locations = '/home/xielab/YiyiWang/FC analysis/Preprocessing/GSN-HydroCel-129.sfp';
    EEG=pop_chanedit(EEG, 'load',{channel_locations 'filetype' 'autodetect'});
    EEG = eeg_checkset( EEG );

% load the eventfile info
ECI_events = load([datapath filename],'ECI*');
fieldname  = fieldnames(ECI_events);
ECI_events = ECI_events.(fieldname{1});
% find the data sample for the baseline section
switch bltype
    case 'eop'
        if age < 60
            ds_start_idx = contains(ECI_events(1,:),{'bas+'});
            ds_end_idx   = contains(ECI_events(1,:),{'bas-'});
        else
            ds_start_idx = contains(ECI_events(1,:),{'eop+'});
            ds_end_idx   = contains(ECI_events(1,:),{'eop-'});
        end
    case 'ecl'
            ds_start_idx = contains(ECI_events(1,:),{'ecl+'});
            ds_end_idx   = contains(ECI_events(1,:),{'ecl-'});
end
ds_start     = ECI_events{4,ds_start_idx};
ds_end       = ECI_events{4,ds_end_idx};

%remove the first and end 2s and extract the baseline section
EEG = pop_select(EEG, 'point',[ds_start+srate*2 ds_end-srate*2]);
% temporarily save the channel locations in a variable
chanlocs = EEG.chanlocs;

%% Preprocessing from here

%% 1. Filtering
% ERPLAB IIR filter
EEG  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff', [1 50], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  8 ); 
% EEGLAB FIR filter

% CleanLine function to suppress the 60 Hz line noise
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:EEG.nbchan] ,'computepower',1,'linefreqs',60,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1)

%% 2. Create an eventlist and do segmentation 
[nch nsamples] = size(EEG.data);  
nseconds = nsamples/srate; 
nepochs  = floor(nseconds/segduration);%1 or 2s
item = 0;bepoch = 0;diff=0; dura=1000*segduration; enable=1;ecode = 1; label = 'TARG';
elistname = strrep(filename,'.mat',['.baseline.erplist_' num2str(segduration) 's.txt']);
if exist([datapath elistname])
    delete([datapath elistname])
end
item = 0;
for i = 1:nepochs;
    onset = segduration*(i-1);
    item  = item +1;
    string=sprintf('%i %i %i %s %f %.2f %.1f 00000000 00000000 %i [    ]',item,bepoch,ecode,label,onset,diff,dura,enable);
    disp(string);
    dlmwrite([datapath elistname],string,'-append','delimiter','');  
end
EEG = pop_importeegeventlist( EEG, [datapath elistname], 'ReplaceEventList', 'on' );
EEG = eeg_checkset( EEG);
% ERPLAB segmentation tool
binlistname = fullfile(folderBase,'Binlists','EmokidsBaselineBinList.txt');
EEG         = pop_binlister( EEG , 'BDF', binlistname, 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG' ); % GUI: 24-Oct-2017 15:26:01
EEG         = pop_epochbin(EEG , [0  dura],  'none');  

%% 3. Find out extraordinarily bad channels using Faster 
% ch129 is the reference
list_properties = channel_properties(EEG, 1:EEG.nbchan, 129);
FASTbadIdx=min_z(list_properties);
FASTbadChans=find(FASTbadIdx==1);
% check if the number of bad channels exceeds the threshold;
badChans = FASTbadChans;
badChans(badChans > 128) = []; %the Cz (channel 129) will be removed later;
if length(badChans) > 18 
    disp(['Participant ' num2str(participantnumber) ' has more than 18 bad channels, so should not be used.']);
    return
end
% Remove the bad channels and Cz;
chan2Bremoved = union(badChans,129);
EEG = pop_select( EEG,'nochannel',chan2Bremoved);

%% 3. artifacts detection and rejection
% 0) ICA pruning
nbchan = 128;
if icapruning
 % Find bad epochs and delete them from dataset
    vol_thrs = [-800 800]; % [lower upper] threshold limit(s) in uV.
    emg_thrs = [-100 30]; % [lower upper] threshold limit(s) in dB.
    emg_freqs_limit = [20 40]; % [lower upper] frequency limit(s) in Hz.
    
    % Find the artifacted epochs across all channels and reject them before doing ICA.
    EEG = pop_eegthresh(EEG,1, 1:EEG.nbchan, vol_thrs(1), vol_thrs(2), EEG.xmin, EEG.xmax,0,0);
    EEG = eeg_checkset(EEG);
    
    % 1         : data type (1: electrode, 0: component)
    % 0         : display with previously marked rejections? (0: no, 1: yes)
    % 0         : reject marked trials? (0: no (but store the  marks), 1:yes)
    
    % Find artifaceted epochs by using power threshold in 20-40Hz frequency band.
    % This method mainly rejects muscle movement (EMG) artifacts.
    EEG = pop_rejspec(EEG, 1,'elecrange', 1:EEG.nbchan, 'method', 'fft', 'threshold', emg_thrs ,'freqlimits', emg_freqs_limit, 'eegplotplotallrej', 0, 'eegplotreject', 0);
    EEG = eeg_checkset(EEG);

    % if there are more than 5 bad channels then reject the epoch.
    EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    badchs4epochs= sum(EEG.reject.rejglobalE,1);
    reject_artifacted_epochs=find(badchs4epochs>4);
    %reject_artifacted_epochs=EEG.reject.rejglobal;
    EEG = pop_rejepoch(EEG, reject_artifacted_epochs, 0);

    % run ica
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'stop', 1E-7, 'interupt','off');
    % EEGICA = EEG;
    % Run adjust to mark bad ICA components
    badICs=[]; EEG_copy =[];
    EEG_copy = EEG;
    % calculate the two EOG channels for SASICA;
    % replace the bad channels in case some of them are EOG channels;
    % however, these channels are not replaced yet for the current EEG_copy;
    nbchan = 128;
    tempdata = zeros(nbchan,EEG_copy.pnts,EEG_copy.trials,'single');
    goodChans = ~ismember([1:nbchan],badChans);
    tempdata(goodChans,:,:) = EEG_copy.data;
    EEG1 = EEG_copy;
    EEG1.data = tempdata;
    EEG1.chanlocs = chanlocs(1:128);
    EEG1.nbchan = nbchan;
    EEG1 = eeg_interp(EEG1,badChans);
    VEOG = mean(EEG1.data([17 14 15 21 22],:,:),1);
    HEOG = mean(EEG1.data([25 32 26 128],:,:),1) - mean(EEG1.data([1 2 8 125],:,:),1);
    
    EEG_copy.data(end+1,:,:) = VEOG;
    EEG_copy.data(end+1,:,:) = HEOG;
    EEG_copy.nbchan = size(EEG_copy.data,1);
    EEG_copy.chanlocs(end+1).labels = 'VEOG';
    EEG_copy.chanlocs(end+1).labels = 'HEOG';
    EEG_copy = eeg_checkset( EEG_copy );
    clear EEG1
    % Define and reject the badICs detected in my grogram using SASICA and Adjust;
    cfg = [];cfg.plotarg = 0;cfg.method = 'adjust';
    [badICs EEG_copy] = DefineRemoveArtificialICAs_EmoKids(cfg,EEG_copy); % if include EEG_copy as an output structure then it will be the one after bad ICs removed;;
    save([SegmentAverageFiles 'Participantnumber ' num2str(participantnumber) '_baseline_badICs.mat'],'badICs');

    % Mark the bad ICs 
    for ic=1:length(badICs)
        EEG.reject.gcompreject(1, badICs(ic))=1;
        EEG = eeg_checkset(EEG);
    end
    ICs2remove=find(EEG.reject.gcompreject)
    EEG = eeg_checkset( EEG );
    EEG = pop_subcomp( EEG, ICs2remove, 0); % remove ICs from dataset
    clear EEG_copy
end

% 1). Find artifaceted epochs by detecting outlier voltage
numChans = EEG.nbchan;
if segduration>1 || age<36
    % maybe use [-150 150] for infants and 2s/3s epochs
    EEG = pop_eegthresh(EEG,1, 1:numChans, -100, 100, EEG.xmin, EEG.xmax,0,0);
else
    EEG = pop_eegthresh(EEG,1, 1:numChans, -100, 100, EEG.xmin, EEG.xmax,0,0);
end
EEG = eeg_checkset( EEG );

% Find artifaceted epochs by using thresholding of frequencies in the data.
EEG = pop_rejspec( EEG, 1,'elecrange', 1:numChans , 'method', 'fft', 'threshold', [-100 30] ,'freqlimits', [20 50], 'eegplotplotallrej', 0, 'eegplotreject', 0);

%% Find the number of artifacted epochs and reject them
EEG = eeg_checkset( EEG );
EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
% if there are more than 5 bad channels then reject the epoch; otherwise, interpolate the bad channels;
% include the first round bad channels;
badchannels = EEG.reject.rejglobalE;
badchs4epochs= sum(badchannels,1);
% max number of bad chs per epoch is 18
maxchan2breplaced = 18 - length(badChans);
reject_artifacted_epochs=find(badchs4epochs>maxchan2breplaced);
badchmatrix = EEG.reject.rejglobalE;
% interpolate for channels; Loop through each epoch, select it, run interp, save data
numEpochs = EEG.trials;
tmpData = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
for e = 1:numEpochs
    % Initialize variables EEGe and EEGe_interp;
    if ismember(e,reject_artifacted_epochs) % this epoch has more than 5 bad channels;
    else
        EEGe = [];
        EEGe_interp = [];
        badChanNum = [];

        %% select only this epoch (e)
        EEGe = pop_selectevent( EEG, 'epoch', e, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
        badChanNum = find(badchmatrix(:,e)); % find which channels are bad for this epoch
        EEGe_interp = eeg_interp(EEGe, badChanNum); % interpolate the bad chans for this epoch
        tmpData(:,:,e) = EEGe_interp.data; % store interpolated data into matrix
        EEG.reject.rejglobal(e)=0;
        EEG.reject.rejglobalE(badChanNum,e) = 0;
    end
end
EEG.data = tmpData;
% reject the bad epochs;
if sum(EEG.reject.rejglobal) == EEG.trials; % all trials are bad and will be removed;
    disp(['Participant ' num2str(participantnumber) ' has all trials bad, so should not be used.']);
    return
else
    EEG = pop_rejepoch( EEG, reject_artifacted_epochs ,0);
end

% 3) Interpolate the overall bad channels that were removed from the beginning;
tempdata = zeros(nbchan,EEG.pnts,EEG.trials,'single');
goodChans = ~ismember([1:nbchan],badChans);
tempdata(goodChans,:,:) = EEG.data;
EEG.data = tempdata;
EEG.chanlocs = chanlocs(1:128);
EEG.nbchan = nbchan;
EEG = eeg_interp(EEG,badChans);

% 4) Re-referencing to average;
EEG = pop_reref(EEG, []);

%% Save the cleaned data in EEGLAB and Fieldtrip format;
% Finaldataset = ['Expt 3 Participantnumber ' num2str(participantnumber) ' Age ' agestring ' Baseline ' bltype '_ica' num2str(icapruning) '.set'];
EEG  = pop_saveset(EEG, 'filename',Finaldataset,'filepath', SegmentAverageFiles,'savemode','onefile');

% %% divide into four bands
% switch bandrange
%     case 'theta'
%         SegmentAverageFiles = fullfile(folderBase,'SegmentAverageFiles',filesep,'theta',filesep);
%     case 'alpha'
%         SegmentAverageFiles = fullfile(folderBase,'SegmentAverageFiles',filesep,'alpha',filesep);
%     case 'beta'
%         SegmentAverageFiles = fullfile(folderBase,'SegmentAverageFiles',filesep,'beta',filesep);
%     case 'gamma'
%         SegmentAverageFiles = fullfile(folderBase,'SegmentAverageFiles',filesep,'gamma',filesep);
%         
% end
