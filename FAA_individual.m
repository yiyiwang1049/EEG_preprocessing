% use fieldtrip and remove trend 1/f, all CoG
function out = FFA_individual(cfg,participantnumber)
if ~isfield(cfg, 'age'),        cfg.age         = 84;    end; age         = cfg.age;
if ~isfield(cfg, 'bltype'),     cfg.bltype      = 'ecl'; end; bltype      = cfg.bltype;
if ~isfield(cfg, 'replacearg'), cfg.replacearg  = 0;     end; replacearg  = cfg.replacearg;
if ~isfield(cfg, 'plotarg'),    cfg.plotarg     = 0;     end; plotarg     = cfg.plotarg;
if ~isfield(cfg, 'segduration'),cfg.segduration = 2;     end; segduration = cfg.segduration;
if ~isfield(cfg, 'icapruning'), cfg.icapruning  = 1;     end; icapruning  = cfg.icapruning;

epochlength = segduration;

global EEG;

folderBase  = '/home/xielab/YiyiWang/FC analysis/Preprocessing';

if cfg.age==36
    SegmentAverageFiles = fullfile(folderBase,'SegmentAverageFiles/3Y',filesep);
end 

if cfg.age==60
    SegmentAverageFiles = fullfile(folderBase,'SegmentAverageFiles/5Y',filesep);
end 

if cfg.age==84
    SegmentAverageFiles = fullfile(folderBase,'SegmentAverageFiles/7Y',filesep);
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

Finaldataset = ['Expt 3 Participantnumber ' num2str(participantnumber) ' Age ' agestring ' Baseline_' bltype '_' num2str(segduration) 's_ica' num2str(icapruning) '.set'];


%load data
EEG = pop_loadset('filename',Finaldataset,'filepath',SegmentAverageFiles);
EEG = eeg_checkset(EEG);

% transform eeg structure to fieldtrip
data = eeglab2fieldtrip(EEG, 'preprocessing', 'none'); % get the target structure
 
%% FFT (% log power)
cfg = [];
cfg.output = "pow"; % return the power-spectra
cfg.channel = data.label; % default = 'all'
cfg.foi = 1:(1/epochlength):30; % set the start point?
cfg.method = "mtmfft"; % multitaper frequency transformation
cfg.taper = "hanning"; % dpss (default) for multitaper

spectra_fft = ft_freqanalysis(cfg,data);

% % channel plot
% chan = 1;
% figure;
% plot(spectra_fft.freq,(spectra_fft.powspctrm(chan,:)),'LineWidth',2);
% title(spectra_fft.label{chan,1});
% xlabel('Frequency (Hz)');
% ylabel('Power (\mu V^2)');
% posterior channel: 66,67,72,77,71,76,84,70,75,83

%% ROI

chF3 = [27,23,19,24,28,20];
chF4 = [4,3,123,124,118,117];


%% remove calculate the average power of frontal alpha
load IAFout.mat
IAF = IAFout(find(IAFout(:,1)==participantnumber),2);

alphamin = round(.8*IAF*2)/2;
alphamax = round(1.2*IAF*2)/2;

freqmin = 2*alphamin-1;
freqmax = 2*alphamax-1;

% calculate the relative pow from the theta band to gamma;
totalpow = sum(spectra_fft.powspctrm(:,5:59),2); % from 3 Hz to the end;
totalbins= length(spectra_fft.freq(:,5:59));
totalpow = repmat(totalpow,[1 totalbins]);
relpow   = 100*(spectra_fft.powspctrm(:,5:59) ./ totalpow);
relfreq  = spectra_fft.freq(:,5:59);

outF3 = mean(mean(spectra_fft.powspctrm(chF3,freqmin:freqmax)));
outF4 = mean(mean(spectra_fft.powspctrm(chF4,freqmin:freqmax)));
out.ln = log(outF4) - log(outF3);
out.ratio = 100*(outF4-outF3)/(outF4+outF3);
out.relative = mean(mean(relpow(chF4,freqmin:freqmax))) - mean(mean(relpow(chF3,freqmin:freqmax)));
out.lnrelative = log(mean(mean(relpow(chF4,freqmin:freqmax)))) - log(mean(mean(relpow(chF3,freqmin:freqmax))));

end






