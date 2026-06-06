% use fieldtrip and remove trend 1/f, all CoG
% posterior area
function out = IAFanalysisDS(cfg,participantnumber)
if ~isfield(cfg, 'age'),        cfg.age         = 84;    end; age         = cfg.age;
if ~isfield(cfg, 'bltype'),     cfg.bltype      = 'eop'; end; bltype      = cfg.bltype;
if ~isfield(cfg, 'replacearg'), cfg.replacearg  = 0;     end; replacearg  = cfg.replacearg;
if ~isfield(cfg, 'plotarg'),    cfg.plotarg     = 0;     end; plotarg     = cfg.plotarg;
if ~isfield(cfg, 'segduration'),cfg.segduration = 2;     end; segduration = cfg.segduration;
if ~isfield(cfg, 'icapruning'), cfg.icapruning  = 1;     end; icapruning  = cfg.icapruning;
if ~isfield(cfg, 'ROI'), cfg.ROI  = "posterior";     end; ROI  = cfg.ROI;

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

if strcmp(ROI,"posterior")
    ch = [54,79,53,61,62,78,86,52,60,67,72,77,85,92,59,66,71,76,84,91,58,65,70,75,83,90,96,64,69,74,82,89,95];
end

%% eye close
Finaldataset = ['Expt 3 Participantnumber ' num2str(participantnumber) ' Age ' agestring ' Baseline_ecl_' num2str(segduration) 's_ica' num2str(icapruning) '.set'];


%load data
EEG = pop_loadset('filename',Finaldataset,'filepath',SegmentAverageFiles);
EEG = eeg_checkset(EEG);

% transform eeg structure to fieldtrip
data = eeglab2fieldtrip(EEG, 'preprocessing', 'none'); % get the target structure
 
% FFT
cfg = [];
cfg.output = "pow"; % return the power-spectra
cfg.channel = data.label; % default = 'all'
cfg.foi = 1:30;
cfg.method = "mtmfft"; % multitaper frequency transformation
cfg.taper = "hanning"; % dpss (default) for multitaper

spectra_fft_ecl = ft_freqanalysis(cfg,data);

%% eye open
Finaldataset = ['Expt 3 Participantnumber ' num2str(participantnumber) ' Age ' agestring ' Baseline_eop_' num2str(segduration) 's_ica' num2str(icapruning) '.set'];


%load data
EEG = pop_loadset('filename',Finaldataset,'filepath',SegmentAverageFiles);
EEG = eeg_checkset(EEG);

% transform eeg structure to fieldtrip
data = eeglab2fieldtrip(EEG, 'preprocessing', 'none'); % get the target structure
 
%% FFT
cfg = [];
cfg.output = "pow"; % return the power-spectra
cfg.channel = data.label; % default = 'all'
cfg.foi = 1:30;
cfg.method = "mtmfft"; % multitaper frequency transformation
cfg.taper = "hanning"; % dpss (default) for multitaper

spectra_fft_eop = ft_freqanalysis(cfg,data);

%% eye close - eye open
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x1-x2';
spectra_fft = ft_math(cfg, spectra_fft_ecl,spectra_fft_eop);


%% find attenuation with remove trend
x=1:30;
% remove the trend: a second way is to create a 1/f model (F) and then remove the model from the data;
y = mean(squeeze(mean(spectra_fft.powspctrm(ch,:),1)),1)';
nfreq = length(y);

freqbins = 1:nfreq;
OneOverF = (1./freqbins); % create 1/f
y_regress = fitlm(OneOverF,y); %regress out 1/f and get resids
y_resids  = y_regress.Residuals.Raw; %pull out raw resids
% figure; plot(x,y,x,y_resids);

%Gaussian fitting:
nfreq_alpha = length([5:12]);
w = gausswin(nfreq_alpha); %nfreq should be the number of frequecy bins in the target frequency band (alpha);
% the length of w should equal to the frequency bins in the band to be fitted;
% alphapow 
alphafreqbins = 5:12;
alphapsd = y_resids(alphafreqbins);
alphapsd_fit = alphapsd.*w;

% all CoG
cog = sum(alphapsd_fit'.* alphafreqbins) / sum(alphafreqbins);
    alphapsd_distance = abs(alphapsd_fit-cog);
    out.CoG = alphafreqbins(find(alphapsd_distance==min(alphapsd_distance)));


[pks0, locs0, widths0, proms0] = findpeaks(alphapsd);
[pks0g, locs0g, widths0g, proms0g] = findpeaks(alphapsd_fit);

out.numpk = length(pks0g);

% pks = y-value for peak
% locs = freq. where peak is located 
% widths = width of peak
% proms = prominance of peak

if length(pks0g) > 1 % two peaks: find the center of gravity
    out.peak1 = alphafreqbins(locs0g(1));
    out.peak2 = alphafreqbins(locs0g(2));
    out.final = out.CoG;
elseif length(pks0g) == 1 % one peak;
    out.peak1 = alphafreqbins(locs0g(1));
    out.peak2 = 0;
    out.final = out.peak1;
else 
    out.peak1 = 0;
    out.peak2 = 0;
    out.final = 0;

end
    out.abspk = max(alphafreqbins(locs0g));
    out.avgpk = mean(alphafreqbins(locs0g));

%% find attenuation with out removing trend
x=1:30;
y = mean(squeeze(mean(spectra_fft.powspctrm(ch,:),1)),1)';
nfreq = length(y);
alphafreqbins = 5:12;
alphapsd = y(alphafreqbins);

[pks0, locs0, widths0, proms0] = findpeaks(alphapsd);

out.attnumpk = length(pks0);
cog = sum(alphapsd'.* alphafreqbins) / sum(alphafreqbins);
    alphapsd_distance = abs(alphapsd-cog);
    out.attCoG = alphafreqbins(find(alphapsd_distance==min(alphapsd_distance)));

% pks = y-value for peak
% locs = freq. where peak is located 
% widths = width of peak
% proms = prominance of peak

if length(pks0) > 1 % two peaks: find the center of gravity
    out.attpeak1 = alphafreqbins(locs0(1));
    out.attpeak2 = alphafreqbins(locs0(2));
    out.attfinal = out.attCoG;
        out.attabspk = max(alphafreqbins(locs0));
    out.attavgpk = mean(alphafreqbins(locs0));
elseif length(pks0) == 1 % one peak;
    out.attpeak1 = alphafreqbins(locs0(1));
    out.attpeak2 = 0;
    out.attfinal = out.attpeak1;
        out.attabspk = max(alphafreqbins(locs0));
    out.attavgpk = mean(alphafreqbins(locs0));
else 
    out.attpeak1 = 0;
    out.attpeak2 = 0;
    out.attfinal = 0;
        out.attabspk = 0;
    out.attavgpk = 0;

end





end






