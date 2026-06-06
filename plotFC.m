% Control
% Default
% Dorsal attention
% Limbic
% Ventral attention
% SomatoMotor
% Visual

load('plotalphaRo.mat');
x = plotRo{1};
x= x';
test.ImageGridAmp = tanh(x);

load('plotthetaRo.mat');
x = plotRo{7};
x= x';
test.ImageGridAmp = tanh(x);

load('plotbetaRo.mat');
x = plotRo{6};
x= x';
test.ImageGridAmp = tanh(x);

%% plot the network
load('network_firstROI_vertices_5yf.mat')
load('networkver_5yf.mat')
load('schaefer2018_5yf.mat'); 

i = 1; networkver_first{i} = schaefer2018_5yf(7).Vertices;
i = 3; networkver_first{i} = schaefer2018_5yf(43).Vertices;
i = 4; networkver_first{i} = schaefer2018_5yf(57).Vertices;
i = 5; networkver_first{i} = schaefer2018_5yf(65).Vertices;

for i =1:7
    position{i} = zeros(5000,1);
    for j=1:5000
        if ismember(j,networkver{i})
            position{i}(j) = 0;
        end
                if ismember(j,networkver_first{i})
            position{i}(j) = 1;
        end

    end
end


test.ImageGridAmp = position{1};

