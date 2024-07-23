% =============================================================================
% Description: Code of the article "Actuation manifold from snapshots data"
% Authors: Luigi Marra, Guy Y. Cornejo Maceda, Andrea Meilán-Vila, Vanesa Guerrero, 
% Salma Rashwan, Bernd R. Noack, Stefano Discetti, and Andrea Ianiro.
% DOI: To be received
% Dataset DOI: To be received
% GitHub: https://github.com/Lmarra1/Actuation-manifold-from-snapshot-data.git
% Description: 
% This script performs ISOMAP dimensionality reduction on a 
% dataset of snapshots of the fluidic pinball with control.
% =============================================================================
%
%
% Dataset and Variables Description:
% -----------------------------------------------------------------------------
% 1. **Grids (X_new and Y_new)**:
%    - X_new and Y_new are meshgrid matrices defining a structured grid. 
%    - The grids have a spacing of dx = 0.05 and dy = 0.05, providing the 
%      coordinates for interpolating snapshots.
%
% 2. **Snapshots**:
%    - Snapshots are interpolated on the structured grid defined by X_new and Y_new.
%    - Snapshots are stored in matrices Snap_u (u-component) and Snap_v (v-component).
%    - The snapshots are compiled into a single matrix U. Each column of U 
%      represents a snapshot with the velocity components u and v concatenated 
%      for every point on the grid.
%
% 3. **Actuation Laws (Matrix b)**:
%    - Matrix b contains the actuation laws used in the simulation. Each row 
%      corresponds to a different cylinder rotation speed: front, top, and bottom.
%    - The actuation laws include all combinations of speeds ranging from -3 to 
%      3 with a step size of 1.
%    - For each actuation law, 20 snapshots are recorded after a transient period 
%      (800 c.u. from the start of the simulation) with 
%      a time step of 1 c.u.
% -----------------------------------------------------------------------------
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
ReadH5(Data_path + "\Grid.h5");                % Load grid data (X_new and Y_new)
ReadH5(Data_path + "\TrainingDataset_1.h5");   % Load training dataset 1 (Additional data)
% ReadH5(Data_path + "\TrainingDataset_2.h5"); % Load training dataset 2 (u-component snapshots)
% ReadH5(Data_path + "\TrainingDataset_3.h5"); % Load training dataset 3 (v-component snapshots)

%% Create snapshots matrix
% Combine u and v components into a single matrix U
U = [Snap_u; Snap_v];
% Apply mask to set velocities at cylinder locations to zero
U(logical([Mask(:); Mask(:)]), :) = 0; 

%% ISOMAP Analysis
% Calculate the distance matrix
D = L2_distance(U, U);          % Compute pairwise L2 distances between snapshots
D = D * sqrt(dA);               % Adjust distances to match discrete definitions
D = (D + D.') / 2;              % Symmetrize the distance matrix
for i = 1:size(D, 1)
    D(i, i) = 0;                % Set diagonal elements to zero
end

%% Perform ISOMAP
ke = 40; % Number of nearest neighbors for ISOMAP (Selected k_e = 40)

ISOoptions.display = 0; % no plot of residual variance and 2-D embedding
ISOoptions.overlay = 0; % no overlay graph on 2-D embedding
ISOoptions.verbose = 1; % display progress report

% Iterate over each k-value (in this case, just one value)
for index = 1:length(ke)

    % Perform ISOMAP dimensionality reduction
    [Y, R, ~, Dg] = IsoMap(D, 'k', ke(index), ISOoptions);
    
    % Store results in structure for further analysis
    ISOMAP.Gamma.("k" + string(ke(index))) = Y.coords{10, 1};   % Low-dimensional embedding coordinates
    ISOMAP.R.("k" + string(ke(index)))     = R;                 % Residual variance of the embedding
    ISOMAP.Dg.("k" + string(ke(index)))    = Dg;                % Geodesic distance matrix
    ISOMAP.Phi.("k" + string(ke(index)))   = U*ISOMAP.Gamma.("k" + string(ke(index)))./(sqrt(sum(Gamma.^2, 2))).'; % Pseudomodes                 
    
    % Display progress for the current k-value
    disp("------------------------>  k: " + string(ke(index)));
end

disp("Calculations completed. You can view the results in the variable ""ISOMAP""")
