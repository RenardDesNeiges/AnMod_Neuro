%% Script for PCA 

addpath(genpath('..'))
clear 
H01_TDM_2kmh_features = load('features/H01_TDM_2kmh_features.mat').s;
H01_TDM_35kmh_features = load('features/H01_TDM_35kmh_features.mat').s;

%% Converting the struct to a vector

H01_TDM_2kmh_vec = cell2mat(struct2cell(H01_TDM_2kmh_features));
H01_TDM_35kmh_vec = cell2mat(struct2cell(H01_TDM_35kmh_features));

X = [H01_TDM_2kmh_vec,H01_TDM_35kmh_vec]';

%% Running PCA on the data

[coeff,score,latent,~,explained] = pca(X);

covarianceMatrix = cov(X);
[V,D] = eig(covarianceMatrix);

imagesc(coeff)