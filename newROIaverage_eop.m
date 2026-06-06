clear;clc;
% select vertices
load('schaefer2018_7yf.mat');
schaefer2018 = schaefer2018_7yf;
for i = 1:100
    selectver{i} = schaefer2018(i+2).Vertices;
end
%% read the 5000*5000
load('thetaRo_eop_newIAF.mat');
thetaRo_eop = Ro;
for num = 1:length(thetaRo_eop)

% read FC values and transform back to correlation
FCvalues = thetaRo_eop{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of Ro_eopIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewthetaRo_eop{num} = FCnew;


end

save("newthetaRo_eop_newIAF.mat","FCnewthetaRo_eop",'-v7.3');

% read the 5000*5000
load('thetaRp_eop_newIAF.mat');
thetaRp_eop = Rp;
for num = 1:length(thetaRp_eop)

% read FC values and transform back to correlation
FCvalues = thetaRp_eop{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of Rp_eopIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewthetaRp_eop{num} = FCnew;


end

save("newthetaRp_eop_newIAF.mat","FCnewthetaRp_eop",'-v7.3');
% 
% %% read the 5000*5000
% load('thetaRo_eop2.mat');
% thetaRo_eop2 = Ro;
% for num = 1:length(thetaRo_eop2)
% 
% % read FC values and transform back to correlation
% FCvalues = thetaRo_eop2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro_eop2Is
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
% FCnewthetaRo_eop2{num} = FCnew;
% num
% 
% 
% end
% 
% save("newthetaRo_eop2.mat","FCnewthetaRo_eop2",'-v7.3');
% 
% % read the 5000*5000
% load('thetaRp_eop2.mat');
% thetaRp_eop2 = Rp;
% for num = 1:length(thetaRp_eop2)% change to 1
% 
% % read FC values and transform back to correlation
% FCvalues = thetaRp_eop2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp_eop2Is
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
% FCnewthetaRp_eop2{num} = FCnew;
% num
% 
% end
% 
% save("newthetaRp_eop2.mat","FCnewthetaRp_eop2",'-v7.3');
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
% load('betaRo_eop.mat');
% betaRo_eop = Ro;
% for num = 1:length(betaRo_eop)
% 
% % read FC values and transform back to correlation
% FCvalues = betaRo_eop{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro_eopIs
% for i = 1:100
%     for j = 1:100
%         target = FCvalues(selectver{i},selectver{j});
%         FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
%     end
% 
% end
% 
% FCnewbetaRo_eop{num} = FCnew;
% 
% 
% end
% 
% save("newbetaRo_eop.mat","FCnewbetaRo_eop",'-v7.3');
% 
% % read the 5000*5000
% load('betaRp_eop.mat');
% betaRp_eop = Rp;
% for num = 1:length(betaRp_eop)
% 
% % read FC values and transform back to correlation
% FCvalues = betaRp_eop{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp_eopIs
% for i = 1:100
%     for j = 1:100
%         target = FCvalues(selectver{i},selectver{j});
%         FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
%     end
% 
% end
% 
% FCnewbetaRp_eop{num} = FCnew;
% 
% 
% end
% 
% save("newbetaRp_eop.mat","FCnewbetaRp_eop",'-v7.3');
% 
% %% read the 5000*5000
% load('betaRo_eop2.mat');
% betaRo_eop2 = Ro;
% for num = 1:length(betaRo_eop2)
% 
% % read FC values and transform back to correlation
% FCvalues = betaRo_eop2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro_eop2Is
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
% FCnewbetaRo_eop2{num} = FCnew;
% num
% 
% 
% end
% 
% save("newbetaRo_eop2.mat","FCnewbetaRo_eop2",'-v7.3');
% 
% % read the 5000*5000
% load('betaRp_eop2.mat');
% betaRp_eop2 = Rp;
% for num = 1:length(betaRp_eop2)% change to 1
% 
% % read FC values and transform back to correlation
% FCvalues = betaRp_eop2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp_eop2Is
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
% FCnewbetaRp_eop2{num} = FCnew;
% num
% 
% end
% 
% save("newbetaRp_eop2.mat","FCnewbetaRp_eop2",'-v7.3');
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
load('alphaRo_eop_newIAF.mat');
alphaRo_eop = Ro;
for num = 1:length(alphaRo_eop)

% read FC values and transform back to correlation
FCvalues = alphaRo_eop{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of Ro_eopIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewalphaRo_eop{num} = FCnew;


end

save("newalphaRo_eop_newIAF.mat","FCnewalphaRo_eop",'-v7.3');

% read the 5000*5000
load('alphaRp_eop_newIAF.mat');
alphaRp_eop = Rp;
for num = 1:length(alphaRp_eop)

% read FC values and transform back to correlation
FCvalues = alphaRp_eop{num};

% calculate the FC matrixs for 100*100 by averaging the correlations for
% every pair of Rp_eopIs
for i = 1:100
    for j = 1:100
        target = FCvalues(selectver{i},selectver{j});
        FCnew(i,j) = nanmean(target(target~=Inf));%remove Inf values
    end

end

FCnewalphaRp_eop{num} = FCnew;


end

save("newalphaRp_eop_newIAF.mat","FCnewalphaRp_eop",'-v7.3');
% 
% %% read the 5000*5000
% load('alphaRo_eop2.mat');
% alphaRo_eop2 = Ro;
% for num = 1:length(alphaRo_eop2)
% 
% % read FC values and transform back to correlation
% FCvalues = alphaRo_eop2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Ro_eop2Is
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
% FCnewalphaRo_eop2{num} = FCnew;
% num
% 
% 
% end
% 
% save("newalphaRo_eop2.mat","FCnewalphaRo_eop2",'-v7.3');
% 
% % read the 5000*5000
% load('alphaRp_eop2.mat');
% alphaRp_eop2 = Rp;
% for num = 1:length(alphaRp_eop2)% change to 1
% 
% % read FC values and transform back to correlation
% FCvalues = alphaRp_eop2{num};
% 
% % calculate the FC matrixs for 100*100 by averaging the correlations for
% % every pair of Rp_eop2Is
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
% FCnewalphaRp_eop2{num} = FCnew;
% num
% 
% end
% 
% save("newalphaRp_eop2.mat","FCnewalphaRp_eop2",'-v7.3');
% 


