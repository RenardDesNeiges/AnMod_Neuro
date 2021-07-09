%% Script for multi-species PCA analysis

% Available datasets : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% | Dataset name                      | Species | Condition               |
% | --------------------------------- | ------- | ------------------------|
% | H01_TDM_2kmh                      | human   | Healthy                 |
% | H01_TDM_35kmh                     | human   | Healthy                 |
% | H01_TDM_2kmh_20_incl              | human   | Healthy                 |
% | DM002_TDM_08_1kmh                 | human   | SCI EES                 |
% | DM002_TDM_08_2kmh                 | human   | SCI EES                 |
% | DM002_TDM_1kmh_NoEES              | human   | SCI No EES              |
% | Elektra_20190425_TM20_004         | NHP     | Healthy                 |
% | Elektra_20190425_TM30_002         | NHP     | Healthy                 |
% | Elektra_20190425_TM40_005         | NHP     | Healthy                 |
% | Healthy_330_BIP_RW_06             | rodent  | Healthy                 |
% | Healthy_330_BIP_RW_07             | rodent  | Healthy                 |
% | SCI_Trained_207_RW_STIM_25_04     | rodent  | SCI EES                 |
% | SCI_Trained_207_RW_STIM_35_02     | rodent  | SCI EES                 |
% | SCI_Trained_207_RW_STIM_BWS40_10  | rodent  | SCI EES                 |
% | SCI_trained_207_RW_SPONT_BWS45_05 | rodent  | SCI No EES              |
% | SCI_trained_207_RW_SPONT_BWS40_03 | rodent  | SCI No EES              |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





addpath(genpath('..'))
clc 
clear 
close all
set(0,'DefaultFigureWindowStyle','docked') 

show_plots = true; %set to true to display results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading all datasets and labeling them
datasets  = [   "H01_TDM_2kmh",...
                "H01_TDM_35kmh",...
                "H01_TDM_2kmh_20_incl",...
                "DM002_TDM_08_1kmh",...
                "DM002_TDM_08_2kmh",...
                "DM002_TDM_1kmh_NoEES",...
                "Elektra_20190425_TM20_004",...
                "Elektra_20190425_TM30_002",...
                "Elektra_20190425_TM40_005",...
                "Healthy_330_BIP_RW_06",...
                "Healthy_330_BIP_RW_07",...
                "SCI_Trained_207_RW_STIM_35_02",...
                "SCI_Trained_207_RW_STIM_BWS40_10",...
                "SCI_trained_207_RW_SPONT_BWS40_03",...
                "SCI_trained_207_RW_SPONT_BWS45_05"];
            
species  = [    "Human",...
                "Human",...
                "Human",...
                "Human",...
                "Human",...
                "Human",...
                "NHP",...
                "NHP",...
                "NHP",...
                "Rodent",...
                "Rodent",...
                "Rodent",...
                "Rodent",...
                "Rodent",...
                "Rodent"];
            
            
condition  = [  "Healthy",...
                "Healthy",...
                "Healthy",...
                "EES",...
                "EES",...
                "NoEES",...
                "Healthy",...
                "Healthy",...
                "Healthy",...
                "Healthy",...
                "Healthy",...
                "EES",...
                "EES",...
                "NoEES",...
                "NoEES"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = [];
for i = 1:1:size(datasets,2)
    name = strcat("features/",datasets(i),"_features.mat");
    struct = load(name).s;
    vec = cell2mat(struct2cell(struct));
    vec(isnan(vec))=0;
    data = [data,vec];
end
% names of the different variables in the feature vector
struct = load(strcat("features/",datasets(1),"_features.mat")).s;
field_names = fieldnames(struct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running PCA on the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zscore(data');

[coefs,score] = pca(X);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Standardized Feature Vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if show_plots == true
    figure
end
for i = 1:1:size(datasets,2)
    vec = X(i,:);
    if show_plots == true
        subplot(5,3,i)
        bar(vec)
         tx = title(datasets(i));
         set(tx,'Interpreter','none')
        set(gcf,'color','w');
    end
end
if show_plots
    sgtitle('Standardized feature vectors for different time series')       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Biplots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if show_plots == true
       %    first 4 eigenvectors biplots
       biplot_dimensions(coefs,score,1:2,condition)
       biplot_dimensions(coefs,score,1:2,species)
       biplot_dimensions(coefs,score,2:3,condition)
       biplot_dimensions(coefs,score,2:3,species)
       biplot_dimensions(coefs,score,3:4,condition)
       biplot_dimensions(coefs,score,3:4,species)

       biplot_dimensions(coefs,score,[3, 1],condition)
       biplot_dimensions(coefs,score,[3, 1],species)

       biplot_dimensions(coefs,score,[4, 1],condition)
       biplot_dimensions(coefs,score,[4, 1],species)

       biplot_dimensions(coefs,score,[4, 2],condition)
       biplot_dimensions(coefs,score,[4, 2],species)
       %    content of the first 4 eigenvectors
       
       bar_eigenvector(coefs,1,field_names)
       bar_eigenvector(coefs,2,field_names)
       bar_eigenvector(coefs,3,field_names)
       bar_eigenvector(coefs,4,field_names)
    
end
