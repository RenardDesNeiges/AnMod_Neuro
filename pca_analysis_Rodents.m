%% Script for PCA analysis of rodent datasets

% loading all datasets

% Available datasets : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% | Dataset name                      | Species | Condition               |
% | --------------------------------- | ------- | ------------------------|
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
datasets  = [   "Healthy_330_BIP_RW_06",...
                "Healthy_330_BIP_RW_07",...
                "SCI_Trained_207_RW_STIM_35_02",...
                "SCI_Trained_207_RW_STIM_BWS40_10",...
                "SCI_trained_207_RW_SPONT_BWS40_03",...
                "SCI_trained_207_RW_SPONT_BWS45_05"];
            
species  = [    "Rodent",...
                "Rodent",...
                "Rodent",...
                "Rodent",...
                "Rodent",...
                "Rodent"];
            
            
condition  = [  "Healthy",...
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
% Biplots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if show_plots == true
       %    first 4 eigenvectors biplots
       biplot_dimensions(coefs,score,1:2,condition)
       biplot_dimensions(coefs,score,2:3,condition)
       biplot_dimensions(coefs,score,3:4,condition)
       biplot_dimensions(coefs,score,[3, 1],condition)
       biplot_dimensions(coefs,score,[4, 1],condition)
       biplot_dimensions(coefs,score,[4, 2],condition)
       
       %    content of the first 4 eigenvectors
       
       bar_eigenvector(coefs,1,field_names)
       bar_eigenvector(coefs,2,field_names)
       bar_eigenvector(coefs,3,field_names)
       bar_eigenvector(coefs,4,field_names)
    
end
