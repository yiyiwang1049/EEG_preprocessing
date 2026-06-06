% use fieldtrip and remove trend 1/f, all CoG
function out = IAFanalysisGch(cfg,participantnumber)
if ~isfield(cfg, 'age'),        cfg.age         = 84;    end; age         = cfg.age;
if ~isfield(cfg, 'bltype'),     cfg.bltype      = 'eop'; end; bltype      = cfg.bltype;
if ~isfield(cfg, 'replacearg'), cfg.replacearg  = 0;     end; replacearg  = cfg.replacearg;
if ~isfield(cfg, 'plotarg'),    cfg.plotarg     = 0;     end; plotarg     = cfg.plotarg;
if ~isfield(cfg, 'segduration'),cfg.segduration = 2;     end; segduration = cfg.segduration;
if ~isfield(cfg, 'icapruning'), cfg.icapruning  = 1;     end; icapruning  = cfg.icapruning;
if ~isfield(cfg, 'ROI'), cfg.ROI  = "posterior";     end; ROI  = cfg.ROI;

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

%% find peaks
if strcmp(ROI,"all")
    ch = 1:128
end
if strcmp(ROI,"TO_L")
    ch = [57,58,63,64,50,59,65,68]
end
if strcmp(ROI,"TO_R")
    ch = [101,95,96,99,100,90,91,94]
end
if strcmp(ROI,"OI")
    ch = [70,71,72,74,75,76,81,82,83]
end
if strcmp(ROI,"CentralZ")
    ch = [7,31,55,80,106]
end
if strcmp(ROI,"Central3")
    ch = [30,36,37,41,42]
end
if strcmp(ROI,"Central4")
    ch = [87,93,103,104,105]
end
if strcmp(ROI,"FrontalZ")
    ch = [5,10,11,12,16,18]
end
if strcmp(ROI,"posterior")
    ch = [54,79,53,61,62,78,86,52,60,67,72,77,85,92,59,66,71,76,84,91,58,65,70,75,83,90,96,64,69,74,82,89,95];
end
% remove the trend: a second way is to create a 1/f model (F) and then remove the model from the data;
x = cfg.foi;%1:(1/epochlength):30;
y = mean(squeeze(mean(spectra_fft.powspctrm(ch,:),1)),1)';
nfreq = length(y);
freqbins = 1:nfreq;
OneOverF = (1./freqbins); % create 1/f
y_regress = fitlm(OneOverF,y); %regress out 1/f and get resids
y_resids  = y_regress.Residuals.Raw; %pull out raw resids
% figure; plot(x,y,x,y_resids);
% figure; plot(x(9:59),y_resids(9:59));
% xlabel("frequency (HZ)");
% ylabel("power (remove 1/f)")

%Gaussian fitting:
nfreq_alpha = length([7:(1/epochlength):12]);
w = gausswin(nfreq_alpha); %nfreq should be the number of frequecy bins in the target frequency band (alpha);
% the length of w should equal to the frequency bins in the band to be fitted;
% alphapow 
if epochlength == 1
    alphafreqbins = 7:12;
    alphapower = 7:12;
else
    allfreq = 1:0.5:30;
    alphafreqbins = [find(allfreq == 7):find(allfreq == 12)];
    alphapower = 7:0.5:12;
   % alphapower = 1:0.5:(1+0.5*(nfreq_alpha-1));
   % alphapower = 1:nfreq_alpha;
end

alphapsd = y_resids(alphafreqbins);
alphapsd_fit = alphapsd.*w;
% figure;plot(7:(1/epochlength):12,alphapsd,7:(1/epochlength):12,alphapsd_fit);

% all CoG (using raw power or remove 1/f?)
out.CoG = sum(y(alphafreqbins)'.* alphapower) / sum(y(alphafreqbins));
% 7.2391(5:12) 7.6873(6:12) 8.2368(7:12)

[pks0, locs0, widths0, proms0] = findpeaks(alphapsd);
[pks0g, locs0g, widths0g, proms0g] = findpeaks(alphapsd_fit);

out.numpk = length(pks0g);

% pks = y-value for peak
% locs = freq. where peak is located 
% widths = width of peak
% proms = prominance of peak

if length(pks0g) > 1 % two peaks: find the center of gravity, compare prominence
    out.peak1 = x(alphafreqbins(locs0g(1)));
    out.peak2 = x(alphafreqbins(locs0g(2)));
    out.abspk = x(max(alphafreqbins(locs0g)));
    out.avgpk = mean(x(alphafreqbins(locs0g)));
    out.prom = alphapower(locs0g(find(proms0g==max(proms0g))));
    out.final = out.prom;
elseif length(pks0g) == 1 % one peak;
    out.peak1 = x(alphafreqbins(locs0g(1)));
    out.peak2 = 0;
    out.abspk = x(max(alphafreqbins(locs0g)));
    out.avgpk = mean(x(alphafreqbins(locs0g)));
    out.prom = out.peak1;
    out.final = out.peak1;
else 
    out.peak1 = 0;
    out.peak2 = 0;
    out.abspk = 0;
    out.avgpk = 0;
    out.prom = 0;
    out.final = out.CoG;

end
    

end






