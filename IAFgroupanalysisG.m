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
num = 0;
IAF = [];

%% individual
for i = 1:length(FileNames)
    if strcmp(FileNames{i,1}(45:47),bltype) || strcmp(FileNames{i,1}(46:48),bltype)|| strcmp(FileNames{i,1}(47:49),bltype)
        num = num+1;
        numend = length(FileNames{i,1})-33
        participantnumber = str2num(FileNames{i,1}(26:numend));% add eyes close or open
        IAF(num,1) = participantnumber;
        IAF(num,2) = IAFanalysisGch(cfg,participantnumber).CoG;
        IAF(num,3) = IAFanalysisGch(cfg,participantnumber).numpk;
        IAF(num,4) = IAFanalysisGch(cfg,participantnumber).peak1;
        IAF(num,5) = IAFanalysisGch(cfg,participantnumber).peak2;
        IAF(num,6) = IAFanalysisGch(cfg,participantnumber).abspk;
        IAF(num,7) = IAFanalysisGch(cfg,participantnumber).avgpk;
        IAF(num,8) = IAFanalysisGch(cfg,participantnumber).prom;
        IAF(num,9) = IAFanalysisGch(cfg,participantnumber).final;
    end
    
end

T = array2table(IAF,'VariableNames',{'id';'CoG';'numpk';'peak1';'peak2';'abspk';'avgpk';'prom';'final'});
writetable(T,[num2str(age) '_' bltype '_' 'IAF.xlsx']);

%% average

%.4*IAF¨C.8*IAF, .8*IAF¨C1.2*IAF, 1.2*IAF-30, and 30¨C50 Hz 
% % for theta, alpha, beta, and gamma, respectively
% theta = [.4*averagevalue,.8*averagevalue]
% alpha = [.8*averagevalue, 1.2*averagevalue]
% beta = [1.2*averagevalue,30]
% gamma = [30,50]

