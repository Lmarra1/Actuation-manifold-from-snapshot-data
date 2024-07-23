% =============================================================================
% Description: Code of the article "Actuation manifold from snapshots data"
% Authors: Luigi Marra, Guy Y. Cornejo Maceda, Andrea Meilán-Vila, Vanesa Guerrero, 
% Salma Rashwan, Bernd R. Noack, Stefano Discetti, and Andrea Ianiro.
% DOI: To be received
% Dataset DOI: To be received
% GitHub: https://github.com/Lmarra1/Actuation-manifold-from-snapshot-data.git
% Description: 
% This script performs the reconstruction of the snapshots in the test
% dataset starting from the predictions of the MLP.
% MLP maps the actuation+sensors information to manifold coordinates. These
% are then used to reconstruct the snapshots using a KNN interpolation.
% =============================================================================


%% Clear previous data and settings
% Clear all figures, variables, and command window
close all;
clear;
clc;
warning off; % Disable warnings


%% Path to data directory
% Update this path to the location where the data downloaded from Zenodo (DOI: [To be received]) is stored
Data_path = "G:\Mi unidad\DRIVE\PhD\Submission JFM ISOMAP\OA dataset";
% Data_path = "SPECIFY YOUR DATA PATH";

% Load grid and snapshot data from HDF5 files
ReadH5(Data_path + "\Grid.h5");                  % Load grid data (X_new and Y_new)
% ReadH5(Data_path + "\TrainingDataset_1.h5");   % Load training dataset 1 (Additional data)
ReadH5(Data_path + "\TrainingDataset_2.h5");     % Load training dataset 2 (u-component snapshots)
ReadH5(Data_path + "\TrainingDataset_3.h5");     % Load training dataset 3 (v-component snapshots)
U = [Snap_u; Snap_v];                            % Snapshot matrix test dataset
ReadH5(Data_path + "\IsomapResults.h5");         % Isomap results for k_e = 40
ReadH5(Data_path + "\TestDataset.h5");           % Load training dataset 3 (v-component snapshots)
U_t = [Snap_u; Snap_v];                          % Snapshot matrix training dataset



%% Initialization
n_coordinates = 5;  % Number of coordinates to use
Npts = numel(X_new);  % Number of points in the new dataset

% Restrict Gamma to the first n_coordinates
Gamma = Gamma(1:n_coordinates,:);

%% Select the points in the latent space that you want to decode
caseMLP = 2;  % Case selection: 1 for sensing information (cl and cl delayed), 2 for sensing info (x component at points (7, 1.25) and (10, 1.25))
k_d = 40;  % Number of nearest neighbors for Isomap

% Initialize U_hat matrix to store the decoded data
U_hat = nan(Npts*2, size(U_t, 2));

% Loop over all 22 cases
for Case = 1:22
    % Read MLP predictions for the current case
    MLP_pred = readmatrix("MLP_predictions/Case_" + string(Case) + "_pred_dataMLP_" + string(caseMLP) + ".txt");
    
    % Extract the predicted gamma values (last n_coordinates columns)
    gamma_pred = MLP_pred(:, end-n_coordinates+1:end).';
    
    % Decode the latent space coordinates to reconstruct U_hat using Isomap
    U_hat(:, (Case-1)*20+1 : Case*20) = isomap_decoder(gamma_pred, U, Gamma, k_d);
    
    % Display progress
    disp("Reconstruction " + string(Case) + "/22 - done")
end

%% Calculate Cosine Similarity
% Initialize the CosSim array
CosSim = nan(1, 440);

% Loop backwards through each set of columns in U_hat and U_t
for i = 440:-1:1
    % Compute cosine similarity between the corresponding columns of U_hat and U_t
    CosSim(i) = dot(U_hat(:, i), U_t(:, i)) / (norm(U_hat(:, i)) * norm(U_t(:, i)));
end










