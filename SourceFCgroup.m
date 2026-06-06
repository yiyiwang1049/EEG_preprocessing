clear;clc;
cfg=[];

% read names
folderBase  = '/home/xielab/YiyiWang/FC analysis/Preprocessing';
folderband = fullfile(folderBase,'SegmentAverageFiles',filesep);
folderdata = fullfile(folderband,'theta',filesep); %theta
                 

File = dir(fullfile(folderdata,'*.set'));  
FileNames = {File.name}';    

%% individual
for i = 1:length(FileNames)
        numend = length(FileNames{i,1})-33;
        participantnumber = str2num(FileNames{i,1}(7:numend));%7:num
        BPCD = SourceAnalysis(cfg,participantnumber); 
        % Use GPU
        [Rp{i}, Ro{i}] = OrthogonalPowCorr(BPCD,0);
end

for i = 1:length(FileNames)
       numend = length(FileNames{i,1})-33;
        participantnumber(i) = str2num(FileNames{i,1}(6:numend));  
end

save("thetaparticipantnumber_ecl_newIAF.mat","participantnumber");
% use participantnumber2

save('thetaRo_ecl_newIAF.mat','Ro','-v7.3');save('thetaRp_ecl_newIAF.mat','Rp','-v7.3');

