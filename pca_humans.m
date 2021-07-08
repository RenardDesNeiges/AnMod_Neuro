%% Script for PCA 

% loading all (human) datasets

% Available datasets : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% | Dataset name         | Species | Condition                            |
% | -------------------- | ------- | -------------------------------------|
% | H01_TDM_2kmh         | Human   | Healthy, 2kmh walk                   |
% | H01_TDM_35kmh        | Human   | Healthy, 3.5kmh walk                 |
% | H01_TDM_2kmh_20_incl | Human   | Healthy, 2kmh walk on a slope (20°)  |
% | DM002_TDM_08_2kmh    | Human   | SCI, EES treatement, 2kmh walk       |
% | DM002_TDM_08_1kmh    | Human   | SCI, EES treatement, 1khm walk       |
% | DM002_TDM_1kmh_NoEES | Human   | SCI, no EES treatement, 1kmh walk    |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('..'))
clc 
clear 
close all
set(0,'DefaultFigureWindowStyle','docked') 

show_plots = true; %set to true to display results

%dataset names
datasets  = [   "H01_TDM_2kmh",...
                "H01_TDM_35kmh",...
                "H01_TDM_2kmh_20_incl",...
                "DM002_TDM_08_2kmh",...
                "DM002_TDM_08_1kmh"];

data = [];
if show_plots == true
    figure
end
for i = 1:1:size(datasets,2)
    name = strcat("features/",datasets(i),"_features.mat");
    struct = load(name).s;
    vec = cell2mat(struct2cell(struct));
    vec(isnan(vec))=0;
    data = [data,vec];
    
    if show_plots == true
        subplot(size(datasets,2),1,i)
        bar(log(vec))
        tx = title(strcat(...
            "log scale feature vector, time series : ",...
            datasets(i))); % avoids interpreting _ as latex
        set(tx,'Interpreter','none')
        set(gcf,'color','w');
    end
end
if show_plots
    sgtitle('Feature vectors for human time series')       
end

% Running PCA on the data

X = zscore(data);


if show_plots == true
    figure
end
for i = 1:1:size(datasets,2)
    vec = X(:,i);
    if show_plots == true
        subplot(size(datasets,2),1,i)
        bar(vec)
        tx = title(strcat(...
            "standardized feature vector, time series : ",...
            datasets(i))); % avoids interpreting _ as latex
        set(tx,'Interpreter','none')
        set(gcf,'color','w');
    end
end
if show_plots
    sgtitle('Standardized feature vectors for human time series')       
end


[coefs,score] = pca(X);
if show_plots == true
    figure
    biplot(coefs(:,1:3),'Scores',score(:,1:3),'VarLabels',datasets);
    set(gcf,'color','w');
end
