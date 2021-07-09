%% Script for human PCA analysis

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
                "DM002_TDM_1kmh_NoEES"];
            
species  = [    "Human",...
                "Human",...
                "Human",...
                "Human",...
                "Human",...
                "Human"];
            
            
condition  = [  "Healthy",...
                "Healthy",...
                "Healthy",...
                "EES",...
                "EES",...
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
