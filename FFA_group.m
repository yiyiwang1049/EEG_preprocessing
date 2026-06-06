clear;clc;
cfg=[];

if ~isfield(cfg, 'age'),        cfg.age         = 84;    end; age         = cfg.age;
if ~isfield(cfg, 'bltype'),     cfg.bltype      = 'ecl'; end; bltype      = cfg.bltype;

% read names
folderBase  = '/home/xielab/YiyiWang/FC analysis/Preprocessing';

if age == 36
SegmentAverageFiles = fullfile(folderBase,'SegmentAverageFiles/3Y',filesep);
end

if age == 60
SegmentAverageFiles = fullfile(folderBase,'SegmentAverageFiles/5Y',filesep);
end

if age == 84
SegmentAverageFiles = fullfile(folderBase,'SegmentAverageFiles/7Y',filesep);
end

File = dir(fullfile(SegmentAverageFiles,'*.set'));  
FileNames = {File.name}';    
num=0;

%% individual FFA
for i = 1:length(FileNames)
     if strcmp(FileNames{i,1}(45:47),bltype) || strcmp(FileNames{i,1}(46:48),bltype)|| strcmp(FileNames{i,1}(47:49),bltype)
        num = num+1;
        numend = length(FileNames{i,1})-33
        participantnumber = str2num(FileNames{i,1}(26:numend));% add eyes close or open
        FFA(num,1) = participantnumber;
        FFA(num,2) = FAA_individual(cfg,participantnumber).ln;
        FFA(num,3) = FAA_individual(cfg,participantnumber).ratio;
        FFA(num,4) = FAA_individual(cfg,participantnumber).relative;
        FFA(num,5) = FAA_individual(cfg,participantnumber).lnrelative;
    end
end

T = array2table(FFA,'VariableNames',{'id';'ln';'ratio';'relative';'lnrelative'});
writetable(T,[num2str(age) '_' bltype '_' 'FAA.xlsx']);
