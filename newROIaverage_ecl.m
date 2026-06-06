clear;clc;
% select vertices
load('schaefer2018_7yf.mat');
schaefer2018 = schaefer2018_7yf;
for i = 1:100
    selectver{i} = schaefer2018(i+2).Vertices;
end
%% read the 5000*5000
load('thetaRo_ecl_newIAF.mat');
thetaRo_ecl = Ro;
for num = 1:length(thetaRo_ecl)

% read FC values and transform back to correlation
FCvalues = thetaRo_ecl{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of Ro_eclIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewthetaRo_ecl{num} = FCnew;


end

save("newthetaRo_ecl_newIAF.mat","FCnewthetaRo_ecl",'-v7.3');

% read the 5000*5000
load('thetaRp_ecl_newIAF.mat');
thetaRp_ecl = Rp;
for num = 1:length(thetaRp_ecl)

% read FC values and transform back to correlation
FCvalues = thetaRp_ecl{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of Rp_eclIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewthetaRp_ecl{num} = FCnew;


end

save("newthetaRp_ecl_newIAF.mat","FCnewthetaRp_ecl",'-v7.3');
% 
% %% read the 5000*5000
% load('thetaRo_ecl2.mat');
% thetaRo_ecl2 = Ro;
% for num = 1:length(thetaRo_ecl2)
% 
% % read FC values and transform back to correlation
% FCvalues = thetaRo_ecl2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro_ecl2Is
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
% FCnewthetaRo_ecl2{num} = FCnew;
% num
% 
% 
% end
% 
% save("newthetaRo_ecl2.mat","FCnewthetaRo_ecl2",'-v7.3');
% 
% % read the 5000*5000
% load('thetaRp_ecl2.mat');
% thetaRp_ecl2 = Rp;
% for num = 1:length(thetaRp_ecl2)% change to 1
% 
% % read FC values and transform back to correlation
% FCvalues = thetaRp_ecl2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp_ecl2Is
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
% FCnewthetaRp_ecl2{num} = FCnew;
% num
% 
% end
% 
% save("newthetaRp_ecl2.mat","FCnewthetaRp_ecl2",'-v7.3');
% 

% %%
% clear;clc;
% % select vertices
% load('schaefer2018_7yf.mat');
% schaefer2018 = schaefer2018_7yf;
% for i = 1:100
%     selectver{i} = schaefer2018(i+2).Vertices;
% end
% 
% %% read the 5000*5000
% load('betaRo_ecl.mat');
% betaRo_ecl = Ro;
% for num = 1:length(betaRo_ecl)
% 
% % read FC values and transform back to correlation
% FCvalues = betaRo_ecl{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro_eclIs
% for i = 1:100
%     for j = 1:100
%         target = FCvalues(selectver{i},selectver{j});
%         FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
%     end
% 
% end
% 
% FCnewbetaRo_ecl{num} = FCnew;
% 
% 
% end
% 
% save("newbetaRo_ecl.mat","FCnewbetaRo_ecl",'-v7.3');
% 
% % read the 5000*5000
% load('betaRp_ecl.mat');
% betaRp_ecl = Rp;
% for num = 1:length(betaRp_ecl)
% 
% % read FC values and transform back to correlation
% FCvalues = betaRp_ecl{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp_eclIs
% for i = 1:100
%     for j = 1:100
%         target = FCvalues(selectver{i},selectver{j});
%         FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
%     end
% 
% end
% 
% FCnewbetaRp_ecl{num} = FCnew;
% 
% 
% end
% 
% save("newbetaRp_ecl.mat","FCnewbetaRp_ecl",'-v7.3');
% 
% %% read the 5000*5000
% load('betaRo_ecl2.mat');
% betaRo_ecl2 = Ro;
% for num = 1:length(betaRo_ecl2)
% 
% % read FC values and transform back to correlation
% FCvalues = betaRo_ecl2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro_ecl2Is
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
% FCnewbetaRo_ecl2{num} = FCnew;
% num
% 
% 
% end
% 
% save("newbetaRo_ecl2.mat","FCnewbetaRo_ecl2",'-v7.3');
% 
% % read the 5000*5000
% load('betaRp_ecl2.mat');
% betaRp_ecl2 = Rp;
% for num = 1:length(betaRp_ecl2)% change to 1
% 
% % read FC values and transform back to correlation
% FCvalues = betaRp_ecl2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp_ecl2Is
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
% FCnewbetaRp_ecl2{num} = FCnew;
% num
% 
% end
% 
% save("newbetaRp_ecl2.mat","FCnewbetaRp_ecl2",'-v7.3');
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
load('alphaRo_ecl_newIAF.mat');
alphaRo_ecl = Ro;
for num = 1:length(alphaRo_ecl)

% read FC values and transform back to correlation
FCvalues = alphaRo_ecl{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of Ro_eclIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewalphaRo_ecl{num} = FCnew;


end

save("newalphaRo_ecl_newIAF.mat","FCnewalphaRo_ecl",'-v7.3');

% read the 5000*5000
load('alphaRp_ecl_newIAF.mat');
alphaRp_ecl = Rp;
for num = 1:length(alphaRp_ecl)

% read FC values and transform back to correlation
FCvalues = alphaRp_ecl{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of Rp_eclIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewalphaRp_ecl{num} = FCnew;


end

save("newalphaRp_ecl_newIAF.mat","FCnewalphaRp_ecl",'-v7.3');
% 
% %% read the 5000*5000
% load('alphaRo_ecl2.mat');
% alphaRo_ecl2 = Ro;
% for num = 1:length(alphaRo_ecl2)
% 
% % read FC values and transform back to correlation
% FCvalues = alphaRo_ecl2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro_ecl2Is
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
% FCnewalphaRo_ecl2{num} = FCnew;
% num
% 
% 
% end
% 
% save("newalphaRo_ecl2.mat","FCnewalphaRo_ecl2",'-v7.3');
% 
% % read the 5000*5000
% load('alphaRp_ecl2.mat');
% alphaRp_ecl2 = Rp;
% for num = 1:length(alphaRp_ecl2)% change to 1
% 
% % read FC values and transform back to correlation
% FCvalues = alphaRp_ecl2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp_ecl2Is
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
% FCnewalphaRp_ecl2{num} = FCnew;
% num
% 
% end
% 
% save("newalphaRp_ecl2.mat","FCnewalphaRp_ecl2",'-v7.3');
% 


