
%%
clear;clc;
% select vertices
load('ROIver.mat');
for i = 1:8
    selectver{i} = ROIver{i};
end

% read the 100***100
load('7Y/ROI/newthetaRo.mat');
thetaRo = FCnewthetaRo;

for num = 1:length(thetaRo)

FCvalues = thetaRo{num};


for i = 1:8
    for j = 1:8
        FCnew(i,j) = nanmean(nanmean(FCvalues(selectver{i},selectver{j})));
    end

end

FCnewthetaRo{num} = FCnew;


end

save("7Y/network/networkthetaRo.mat","FCnewthetaRo",'-v7.3');

% read the 100***100
load('7Y/ROI/newthetaRp.mat');
thetaRp = FCnewthetaRp;

for num = 1:length(thetaRp)

FCvalues = thetaRp{num};


for i = 1:8
    for j = 1:8
        FCnew(i,j) = nanmean(nanmean(FCvalues(selectver{i},selectver{j})));
    end

end

FCnewthetaRp{num} = FCnew;


end

save("7Y/network/networkthetaRp.mat","FCnewthetaRp",'-v7.3');


%%
clear;clc;
% select vertices
load('ROIver.mat');
for i = 1:8
    selectver{i} = ROIver{i};
end

% read the 100***100
load('7Y/ROI/newalphaRo.mat');
alphaRo = FCnewalphaRo;

for num = 1:length(alphaRo)

FCvalues = alphaRo{num};


for i = 1:8
    for j = 1:8
        FCnew(i,j) = nanmean(nanmean(FCvalues(selectver{i},selectver{j})));
    end

end

FCnewalphaRo{num} = FCnew;


end

save("7Y/network/networkalphaRo.mat","FCnewalphaRo",'-v7.3');

% read the 100***100
load('7Y/ROI/newalphaRp.mat');
alphaRp = FCnewalphaRp;

for num = 1:length(alphaRp)

FCvalues = alphaRp{num};


for i = 1:8
    for j = 1:8
        FCnew(i,j) = nanmean(nanmean(FCvalues(selectver{i},selectver{j})));
    end

end

FCnewalphaRp{num} = FCnew;


end

save("7Y/network/networkalphaRp.mat","FCnewalphaRp",'-v7.3');



%%
clear;clc;
% select vertices
load('ROIver.mat');
for i = 1:8
    selectver{i} = ROIver{i};
end

% read the 100***100
load('7Y/ROI/newbetaRo.mat');
betaRo = FCnewbetaRo;

for num = 1:length(betaRo)

FCvalues = betaRo{num};


for i = 1:8
    for j = 1:8
        FCnew(i,j) = nanmean(nanmean(FCvalues(selectver{i},selectver{j})));
    end

end

FCnewbetaRo{num} = FCnew;


end

save("7Y/network/networkbetaRo.mat","FCnewbetaRo",'-v7.3');

% read the 100***100
load('7Y/ROI/newbetaRp.mat');
betaRp = FCnewbetaRp;

for num = 1:length(betaRp)

FCvalues = betaRp{num};


for i = 1:8
    for j = 1:8
        FCnew(i,j) = nanmean(nanmean(FCvalues(selectver{i},selectver{j})));
    end

end

FCnewbetaRp{num} = FCnew;


end

save("7Y/network/networkbetaRp.mat","FCnewbetaRp",'-v7.3');



