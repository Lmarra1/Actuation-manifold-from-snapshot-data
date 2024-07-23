% =============================================================================
% Description: Code of the article "Actuation manifold from snapshots data"
% Authors: Luigi Marra, Guy Y. Cornejo Maceda, Andrea Meilán-Vila, Vanesa Guerrero, 
% Salma Rashwan, Bernd R. Noack, Stefano Discetti, and Andrea Ianiro.
% DOI: To be received
% Dataset DOI: To be received
% GitHub: https://github.com/Lmarra1/Actuation-manifold-from-snapshot-data.git
% Description: 
% This script performs the reconstruction of the snapshots in the test
% dataset using the two step KNN regression.
% KNN is used to obtain from actuation+sensors information to manifold coordinates. These
% are then used to reconstruict the snapshots using a KNN interpolation.
% First KNN step: actuation+sensors information --> manifold coordinates
% Second KNN step: manifold coordinates --> snapshots
% =============================================================================


%% Clear previous data and settings
% Clear all figures, variables, and command window
close all;
clear;
clc;
warning off; % Disable warnings

%% Set default figure properties
fs         = 10;                                     % Font size for plots
width      = 15;                                     % Figure width
PointSize  = 35;                                     % Size of scatter plots points
set(0, 'DefaultFigureColor', 'w');                   % White figure background
set(0, 'DefaultLineLineWidth', 1);                   % Default line width
set(0, 'DefaultAxesFontSize', fs);                   % Font size for axes
set(0, 'DefaultLegendFontSize', fs);                 % Font size for legend
set(0, 'DefaultTextFontSize', fs);                   % Font size for text
set(0, 'DefaultAxesLineWidth', 1);                   % Line width for axes
set(0, 'DefaultAxesBox', 'on');                      % Box around axes
set(0, 'DefaultAxesColor', [1 1 1]);                 % Axes background color
set(0, 'DefaultLegendInterpreter', 'latex');         % LaTeX interpreter for legend
set(0, 'DefaultTextInterpreter', 'latex');           % LaTeX interpreter for text
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');  % LaTeX interpreter for tick labels


%% Path to data directory
% Update this path to the location where the data downloaded from Zenodo (DOI: [To be received]) is stored
Data_path = "G:\Mi unidad\DRIVE\PhD\Submission JFM ISOMAP\OA dataset";
% Data_path = "SPECIFY YOUR DATA PATH";

% Load grid and snapshot data from HDF5 files
ReadH5(Data_path + "\TestDataset.h5");         % Load training dataset 3 (v-component snapshots)
p_t          = p;                % Kiki parameters test dataset
b_t          = b;                % Actuation vector test dataset
Cl_t         = Cl;               % Cl test dataset
Cl_delayed_t = Cl_delayed;       % Dleayed l for test dataset
v7_125_t     = v7_125;           % v(8,1.25) for test dataset
v10_125_t    = v10_125;          % v(10,1.25) for test dataset
U_t          = [Snap_u; Snap_v]; % Snapshots matrix test dataset

ReadH5(Data_path + "\TrainingDataset_1.h5");   % Load training dataset 1 (Additional data)
ReadH5(Data_path + "\Grid.h5");                % Load grid data (X_new and Y_new)
ReadH5(Data_path + "\IsomapResults.h5");       % Isomap results for k_e = 40
ReadH5(Data_path + "\TrainingDataset_2.h5");   % Load training dataset 2 (u-component snapshots)
ReadH5(Data_path + "\TrainingDataset_3.h5");   % Load training dataset 3 (v-component snapshots)
U = [Snap_u; Snap_v];                          % Snapshots matrix training dataset


nPts  = size(U,1)/2;
nSnap = size(U,2);

%% Extract sensor information
CaseSensor = 1;

if CaseSensor == 1
    s1 = Cl;
    s2 = Cl_delayed;
elseif CaseSensor == 2
    s1 = v7_125;
    s2 = v10_125;
end


%% Use LOOCV to select Kd1 and Kd2
clc; close all;

kd1v = 1:100;
ErrGamma  = zeros(1,length(kd1v));
Input_m  = [p; s1; s2];
Output_m = Gamma(1:5,:);


%%%--------------------- select kd1 ------ KNN1:  act+sens -> Manifold coord
for i = 1:length(kd1v)
    if exist("Dist"); clear Dist; end
    indStor = 0;
    kd = kd1v(i);
    
    for j = 1:343
        ActInd = (j-1)*20+1:j*20;
        
        for k = 1:20
            
            Inputi   = Input_m(:,ActInd(k)); Outputi = Output_m(:,ActInd(k));
            
            if ActInd(1) == 1
                Indexes = ActInd(end)+1:nSnap;
                Dist = sqrt(sum((Inputi - Input_m(:,Indexes)).^2,1));
                
            elseif ActInd(end) == nSnap
                Indexes = 1:ActInd(1)-1;
                Dist = sqrt(sum((Inputi - Input_m(:,Indexes)).^2,1));
                
            else
                Indexes = [1:ActInd(1)-1, ActInd(end)+1:nSnap];
                Dist = sqrt(sum((Inputi - Input_m(:,Indexes)).^2,1));
            end
            
            
            [SortedDist, indSort] = sort(Dist);
            IndKNN      = Indexes(indSort(1:kd));
            SortedDist  = SortedDist(1:kd);  
            KnnOutput   = Output_m(:,IndKNN);         
            PredOutput  = sum(KnnOutput./SortedDist,2)/sum(1./SortedDist);
            ErrGamma(i) = ErrGamma(i) + sum( (PredOutput - Outputi).^2 );
            
        end
    end
    disp(string(i) + "/" + length(kd1v))
end




%%%--------------------- select kd2   ------ KNN2: Manifold coord -> Snapshots
kd2v = 1:100;
ErrSnap  = zeros(1,length(kd2v));
Output_m = U;
Input_m  = Gamma(1:5,:);

for i = 1:length(kd2v)
    if exist("Dist"); clear Dist; end
    indStor = 0;
    kd = kd2v(i);
    
    
    for j = 1:343
        
        ActInd = (j-1)*20+1:j*20;
        
        
        for k = 1:20
            
            Inputi   = Input_m(:,ActInd(k)); Outputi = Output_m(:,ActInd(k));
            
            if ActInd(1) == 1
                Indexes = ActInd(end)+1:nSnap;
                Dist = sqrt(sum((Inputi - Input_m(:,Indexes)).^2,1));
                
            elseif ActInd(end) == nSnap
                Indexes = 1:ActInd(1)-1;
                Dist = sqrt(sum((Inputi - Input_m(:,Indexes)).^2,1));
                
            else
                Indexes = [1:ActInd(1)-1, ActInd(end)+1:nSnap];
                Dist = sqrt(sum((Inputi - Input_m(:,Indexes)).^2,1));
            end
            
            
            [SortedDist, indSort] = sort(Dist);    
            IndKNN = Indexes(indSort(1:kd));
            SortedDist = SortedDist(1:kd);
            KnnOutput = Output_m(:,IndKNN);
            PredOutput   = sum(KnnOutput./SortedDist,2)/sum(1./SortedDist);
            ErrSnap(i) = ErrSnap(i) + sum( (PredOutput - Outputi).^2 );
        end
        
    end
    disp(string(i) + "/" + length(kd1v))
end


%% kd1 and kd2 selection
[~,kd1] = min(ErrGamma);
kd1 = kd1v(kd1);

[~,kd2] = min(ErrSnap);
kd2 = kd2v(kd2);

%% Plot of LOOCV errors
close all
ind = 0;

ind = ind+1;
figure(ind)
plot(kd1v,ErrGamma)
xlabel("$kd_1$")
ylabel("Estimation error LOOCV")
title("Actuation + sensors --\> Low-dim embeddings")
xline(kd1,"--")

ind = ind+1;
figure(ind)
plot(kd2v,ErrSnap)
xlabel("$kd_2$")
ylabel("Estimation error LOOCV")
title("Low-dim embeddings --\> Snapshots")
xline(kd2,"--")


%% Estimate the test cases with knn
if CaseSensor == 2
    s1_t = v7_125_t;
    s2_t = v10_125_t;
else
    s1_t = Cl_t;
    s2_t = Cl_delayed_t;
end


% Prepare vars
ncases = size(U_t,2);
Gamma_tpred = nan(5,ncases);
U_tpred = nan(size(U_t,1),ncases);


% First step of decoding
Input_t = [p_t; s1_t; s2_t];
Input_m  = [p; s1; s2];
Output_m = Gamma(1:5,:);

for i=1:ncases
    Inputi = Input_t(:,i);
    Dist = sqrt(sum((Inputi - Input_m).^2,1));
    [SortedDist, indSort] = sort(Dist);
    PredOutput   = sum(Output_m(:,indSort(1:kd1))./SortedDist(1:kd1),2)/sum(1./SortedDist(1:kd1));
    Gamma_tpred(:,i) = PredOutput;
    disp(string(i) + "/" + string(ncases) + " predicted")
end


% Second step of decoding
Input_t  = Gamma_tpred;
Input_m  = Gamma(1:5,:);
Output_m = U;

for i=1:ncases
    Inputi = Input_t(:,i);
    Dist = sqrt(sum((Inputi - Input_m).^2,1));
    [SortedDist, indSort] = sort(Dist);
    PredOutput   = sum(Output_m(:,indSort(1:kd2))./SortedDist(1:kd2),2)/sum(1./SortedDist(1:kd2));
    if max(SortedDist(1:kd2)) <1e-9
       PredOutput =  Output_m(:,indSort(1));
    end
    
    if nnz(isnan(PredOutput)) ~=0
        keyboard
    end
    U_tpred(:,i) = PredOutput;
    disp(string(i) + "/" + string(ncases) + " predicted")
end


%% Calculate errors
CosSim = nan(1,size(U_tpred,2));
for i = 1:size(U_tpred,2)
CosSim(i) = dot(U_tpred(:,i), U_t(:,i))/(norm(U_tpred(:,i))*norm(U_t(:,i)));
end

disp("-----> Predicted snapshots in variable U_tpred")

