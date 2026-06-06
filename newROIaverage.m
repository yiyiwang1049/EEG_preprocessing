clear;clc;
% select vertices
load('schaefer2018_7yf.mat');
schaefer2018 = schaefer2018_7yf;
for i = 1:100
    selectver{i} = schaefer2018(i+2).Vertices;
end

%% read the 5000*5000
load('thetaRo.mat');
thetaRo = Ro;
for num = 1:length(thetaRo)

% read FC values and transform back to correlation
FCvalues = thetaRo{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of RoIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewthetaRo{num} = FCnew;


end

save("newthetaRo.mat","FCnewthetaRo",'-v7.3');

% read the 5000*5000
load('thetaRp.mat');
thetaRp = Rp;
for num = 1:length(thetaRp)

% read FC values and transform back to correlation
FCvalues = thetaRp{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of RpIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewthetaRp{num} = FCnew;


end

save("newthetaRp.mat","FCnewthetaRp",'-v7.3');
% 
% %% read the 5000*5000
% load('thetaRo2.mat');
% thetaRo2 = Ro;
% for num = 1:length(thetaRo2)
% 
% % read FC values and transform back to correlation
% FCvalues = thetaRo2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro2Is
% for i = 1:100
%     for j = 1:100
%         if isnan(FCvalues)
%         FCnew(i,j) = NaN;
%         else
%         target = FCvalues(selectver{i},selectver{j});
%         FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf value
%         end
%     end
% 
% end
% 
% FCnewthetaRo2{num} = FCnew;
% num
% 
% 
% end
% 
% save("newthetaRo2.mat","FCnewthetaRo2",'-v7.3');
% 
% % read the 5000*5000
% load('thetaRp2.mat');
% thetaRp2 = Rp;
% for num = 1:length(thetaRp2)% change to 1
% 
% % read FC values and transform back to correlation
% FCvalues = thetaRp2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp2Is
% for i = 1:100
%     for j = 1:100
%          if FCvalues==0
%          FCnew(i,j) = NaN;
%          else
%          target = FCvalues(selectver{i},selectver{j});
%          FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
%          end
%      end
% 
% end
% 
% FCnewthetaRp2{num} = FCnew;
% num
% 
% end
% 
% save("newthetaRp2.mat","FCnewthetaRp2",'-v7.3');
% 

%%
clear;clc;
% select vertices
load('schaefer2018_7yf.mat');
schaefer2018 = schaefer2018_7yf;
for i = 1:100
    selectver{i} = schaefer2018(i+2).Vertices;
end

%% read the 5000*5000
load('betaRo.mat');
betaRo = Ro;
for num = 1:length(betaRo)

% read FC values and transform back to correlation
FCvalues = betaRo{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of RoIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewbetaRo{num} = FCnew;


end

save("newbetaRo.mat","FCnewbetaRo",'-v7.3');

% read the 5000*5000
load('betaRp.mat');
betaRp = Rp;
for num = 1:length(betaRp)

% read FC values and transform back to correlation
FCvalues = betaRp{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of RpIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewbetaRp{num} = FCnew;


end

save("newbetaRp.mat","FCnewbetaRp",'-v7.3');
% 
% %% read the 5000*5000
% load('betaRo2.mat');
% betaRo2 = Ro;
% for num = 1:length(betaRo2)
% 
% % read FC values and transform back to correlation
% FCvalues = betaRo2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro2Is
% for i = 1:100
%     for j = 1:100
%         if isnan(FCvalues)
%         FCnew(i,j) = NaN;
%         else
%         target = FCvalues(selectver{i},selectver{j});
%         FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf value
%         end
%     end
% 
% end
% 
% FCnewbetaRo2{num} = FCnew;
% num
% 
% 
% end
% 
% save("newbetaRo2.mat","FCnewbetaRo2",'-v7.3');
% 
% % read the 5000*5000
% load('betaRp2.mat');
% betaRp2 = Rp;
% for num = 1:length(betaRp2)% change to 1
% 
% % read FC values and transform back to correlation
% FCvalues = betaRp2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp2Is
% for i = 1:100
%     for j = 1:100
%          if FCvalues==0
%          FCnew(i,j) = NaN;
%          else
%          target = FCvalues(selectver{i},selectver{j});
%          FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
%          end
%      end
% 
% end
% 
% FCnewbetaRp2{num} = FCnew;
% num
% 
% end
% 
% save("newbetaRp2.mat","FCnewbetaRp2",'-v7.3');
% 

%%
clear;clc;
% select vertices
load('schaefer2018_7yf.mat');
schaefer2018 = schaefer2018_7yf;
for i = 1:100
    selectver{i} = schaefer2018(i+2).Vertices;
end

%% read the 5000*5000
load('alphaRo.mat');
alphaRo = Ro;
for num = 1:length(alphaRo)

% read FC values and transform back to correlation
FCvalues = alphaRo{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of RoIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewalphaRo{num} = FCnew;


end

save("newalphaRo.mat","FCnewalphaRo",'-v7.3');

% read the 5000*5000
load('alphaRp.mat');
alphaRp = Rp;
for num = 1:length(alphaRp)

% read FC values and transform back to correlation
FCvalues = alphaRp{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of RpIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewalphaRp{num} = FCnew;


end

save("newalphaRp.mat","FCnewalphaRp",'-v7.3');
% 
% %% read the 5000*5000
% load('alphaRo2.mat');
% alphaRo2 = Ro;
% for num = 1:length(alphaRo2)
% 
% % read FC values and transform back to correlation
% FCvalues = alphaRo2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro2Is
% for i = 1:100
%     for j = 1:100
%         if isnan(FCvalues)
%         FCnew(i,j) = NaN;
%         else
%         target = FCvalues(selectver{i},selectver{j});
%         FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf value
%         end
%     end
% 
% end
% 
% FCnewalphaRo2{num} = FCnew;
% num
% 
% 
% end
% 
% save("newalphaRo2.mat","FCnewalphaRo2",'-v7.3');
% 
% % read the 5000*5000
% load('alphaRp2.mat');
% alphaRp2 = Rp;
% for num = 1:length(alphaRp2)% change to 1
% 
% % read FC values and transform back to correlation
% FCvalues = alphaRp2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp2Is
% for i = 1:100
%     for j = 1:100
%          if FCvalues==0
%          FCnew(i,j) = NaN;
%          else
%          target = FCvalues(selectver{i},selectver{j});
%          FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
%          end
%      end
% 
% end
% 
% FCnewalphaRp2{num} = FCnew;
% num
% 
% end
% 
% save("newalphaRp2.mat","FCnewalphaRp2",'-v7.3');
% 
clear;clc;
% select vertices
load('schaefer2018_7yf.mat');
schaefer2018 = schaefer2018_7yf;
for i = 1:100
    selectver{i} = schaefer2018(i+2).Vertices;
end

