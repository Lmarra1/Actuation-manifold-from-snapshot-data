% =============================================================================
% Description: Code of the article "Actuation manifold from snapshots data"
% Authors: Luigi Marra, Guy Y. Cornejo Maceda, Andrea Meilán-Vila, Vanesa Guerrero,
% Salma Rashwan, Bernd R. Noack, Stefano Discetti, and Andrea Ianiro.
% DOI: To be received
% Dataset DOI: To be received
% GitHub: https://github.com/Lmarra1/Actuation-manifold-from-snapshot-data.git
% Description:
% This script performs the plot of the results of ISOMAP dimensionality reduction on a
% dataset of snapshots of the fluidic pinball with control. Plots of the
% reconstruction
% =============================================================================
%
%
% -----------------------------------------------------------------------------
% This code is used to plot the results of the reconstruction of snapshots 
% in the test dataset. The reconstruction is performed in two different steps 
% and in two different ways. The case indicated with 'MLP' involves 
% reconstruction starting from the vector of actuation parameters and sensor 
% information to coordinates in the low-dimensional manifold using a 
% multi-layer perceptron. Subsequently, the coordinates in the manifold 
% are used to reconstruct the snapshots via KNN interpolation.
%
% The second methodology, indicated with 'KNN', uses KNN regression with 
% distance-weighted averages in both of the aforementioned steps.
%
% Case 1 refers to reconstruction using sensor information from Cl and 
% Cl delayed by a quarter of the average shedding period. 
% Case 2, on the other hand, utilizes the v velocity component at the positions 
% (7, 1.25) and (10, 1.25).
% =============================================================================


%% Clear previous data and settings
close all; clear; clc; warning off; % Disable warnings

%% Set default figure properties
fs         = 10;                                     % Font size for plots
width      = 15;                                     % Figure width
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


%% Load data for plotting
% Update this path to the location where the data downloaded from Zenodo (DOI: [To be received]) is stored
Data_path = "G:\Mi unidad\DRIVE\PhD\Submission JFM ISOMAP\OA dataset";
% Data_path = "H:\Mi unidad\DRIVE\PhD\Submission JFM ISOMAP\OA dataset";
% Data_path = "SPECIFY YOUR DATA PATH";

Fig_path = "ArticleFigures";  % Directory for saving figures
if ~exist(Fig_path, 'dir')
    mkdir(Fig_path)
end

% Check if the toolbox for LIC is available
if ~exist("toolbox_image",'dir')
    disp("To plot the LIC of the pseudomodes you need to download the toolbox_image:")
    disp("Gabriel Peyre (2024). Toolbox image" +...
        "(https://www.mathworks.com/matlabcentral/fileexchange/16201-toolbox-image),"+...
        "MATLAB Central File Exchange.")
else
    addpath(genpath("toolbox_image"));
end

% Load all data
ReadH5(Data_path + "\Grid.h5");                       % Load the grid data
ReadH5(Data_path + "\TrainingDataset_1.h5");          % Load the training dataset
b_train = b(:,20:20:end);                             % Actuation inputs for training dataset
ReadH5(Data_path + "\TestDataset.h5");                % Load the test dataset

% Select the decoding method
DecodingMethod = "MLP";                               % Decide between "KNN" and "MLP"
if ~strcmp(DecodingMethod, "KNN") && ~strcmp(DecodingMethod, "MLP")
    DecodingMethod = "MLP";
    disp("Warning: MLP decoding method selected. Choose between MLP and KNN");
end

% Load reconstructed snapshots
% case 1
ReadH5(Data_path + "\Reconstruction_" + DecodingMethod + "1.h5");        % Reconstruction case 1
U_hat_1 = [Snap_u_hat; Snap_v_hat];

% case 2
ReadH5(Data_path + "\Reconstruction_" + DecodingMethod + "2.h5");        % Reconstruction case 2
U_hat_2 = [Snap_u_hat; Snap_v_hat];

% true snapshots
U_t = [Snap_u; Snap_v];             % Matrix of true snapshots
b_test = b(:,20:20:end);            % Actuation inputs for test dataset

nSnap = size(U_t,2);                % Number of snapshots in the test case
nCase = nSnap/20;                   % Number of actuation cases explored (20 snapshots sampled in post transient for each case)


%% Other general variables
% Define circumference and coordinates for plotting the fluidic pinball cylinders
theta = 2*pi*(0:100)/100;           % Define theta for circumference
xCyl1 = -1.299 + 0.5*cos(theta);    % X coordinate for cylinder 1
yCyl1 = 0.5*sin(theta);             % Y coordinate for cylinder 1
xCyl  = 0.5*cos(theta);             % X coordinate for cylinder 2 and 3
yCyl2 = 0.75 + 0.5*sin(theta);      % Y coordinate for cylinder 2
yCyl3 = -0.75 + 0.5*sin(theta);     % Y coordinate for cylinder 3



% Compute cosine similarity between true and reconstructed snapshts
for i = nSnap:-1:1
    CosSim_1(i) = dot(U_hat_1(:,i),U_t(:,i))./(norm(U_hat_1(:,i))*norm(U_t(:,i)));
    CosSim_2(i) = dot(U_hat_2(:,i),U_t(:,i))./(norm(U_hat_2(:,i))*norm(U_t(:,i)));
end

% Calculate distance to the nearest actuation vector in the training dataset
for i = nCase:-1:1
    nn = sqrt(sum((b_train - b_test(:,i)).^2, 1));
    dist(i) = min(nn);
end
dist = reshape(repmat(dist(:),1,20), 1,[]);
[~,inds] = sort(dist);


% Select indices for visualization
transp = 0.5;
if strcmp(DecodingMethod, "KNN"); ind11  = 276; ind12 = 266; else;  ind11  = 83; ind12 = 200; end
[~, ind21] = min(CosSim_1);
[~, ind22] = min(CosSim_2);

% Generate Line Integral Convolution (LIC) for visualization
xsq = linspace(-5,20,300);
ysq = linspace(-5,5 ,300);
[Xsq,Ysq] = meshgrid(xsq,ysq);

for i = 1:8
    if i == 1
        Ufield = U_t(1:125561,ind11); Vfield = U_t(125562:end,ind11);
    elseif i == 2
        Ufield = U_hat_1(1:125561,ind11); Vfield = U_hat_1(125562:end,ind11);
    elseif i == 3
        Ufield = U_t(1:125561,ind21); Vfield = U_t(125562:end,ind21);
    elseif i == 4
        Ufield = U_hat_1(1:125561,ind21); Vfield = U_hat_1(125562:end,ind21);
    elseif i == 5
        Ufield = U_t(1:125561,ind12); Vfield = U_t(125562:end,ind12);
    elseif i == 6
        Ufield = U_hat_2(1:125561,ind12); Vfield = U_hat_2(125562:end,ind12);
    elseif i == 7
        Ufield = U_t(1:125561,ind22); Vfield = U_t(125562:end,ind22);
    elseif i == 8
        Ufield = U_hat_2(1:125561,ind22); Vfield = U_hat_2(125562:end,ind22);
    end
    Mask = Mask(:);
    Ufield(logical(Mask(1:125561))) = 1; Vfield(logical(Mask(1:125561))) = 1;
    Ut = reshape(Ufield,size(X_new,1),size(X_new,2));
    Vt = reshape(Vfield,size(X_new,1),size(X_new,2));

    Ut = interp2(X_new,Y_new,Ut,Xsq,Ysq,"linear");
    Vt = interp2(X_new,Y_new,Vt,Xsq,Ysq,"linear");

    VelF = cat(3, Vt, Ut);

    VelF = VelF./sqrt(sum(VelF.^2,3));

    % Length of the convolution
    w = 80;

    % Set options
    options = struct();
    options.spot_size = 5;        % Set the size of the features

    % Perform Line Integral Convolution
    M(:,:,i) = perform_lic(VelF, w, options); % Matrix with lic plot
    disp("LIC snapshot " + string(i) + " performed")
end


%% Plotting
close all

saveFig = false;

hf = figure('units','centimeters','position',[10, 3, width, 4.5]);

% Case 1
axes("Position",[0.08 0.26 0.2 0.55]);
for i = 22:-1:1
    ind = (i-1)*20+1:i*20;
    plot(dist(i),CosSim_1(ind),'k.'); hold on
    MeanCorr2(i) = mean(CosSim_1(ind));
end
plot(dist(floor((ind11-1)/20)+1),CosSim_1(ind11),'o',"color",[0.4660 0.6740 0.1880]);
plot(dist(floor((ind21-1)/20)+1),CosSim_1(ind21),'sq',"color",[0.4660 0.6740 0.1880]);
xlim([0 0.8]); ylim([0.8, 1.002]) 
xlabel('$|| \textbf{\emph b}_i - \textbf{\emph b}_{i(1)} ||_2$', 'Interpreter', 'latex');
ylabel("$S_C(\textbf{\emph u}_i,\hat{\textbf{\emph u}}_i)$")
set(gca,"XTick",0:0.2:0.8); set(gca,"XTickLabel",string(0:0.2:0.8))
set(gca,"YTick",0:0.05:1);  set(gca,"YTickLabel",string(0:0.05:1)); 
grid on


% Case 2
axes("Position",[0.725 0.26 0.2 0.55]);
for i = 22:-1:1
    ind = (i-1)*20+1:i*20;
    plot(dist(i),CosSim_2(ind),'k.');hold on
    MeanCorr2(i) = mean(CosSim_2(ind));
end

plot(dist(floor((ind12-1)/20)+1),CosSim_2(ind12),'o',"color",[0.4660 0.6740 0.1880]);
plot(dist(floor((ind22-1)/20)+1),CosSim_2(ind22),'sq',"color",[0.4660 0.6740 0.1880]);
annotation('textbox', [0.1, 0.81 0.1, 0.1], 'String', 'median', 'Color', [0 0 0],"EdgeColor","none","fontsize",9,"Interpreter","latex");
annotation('textbox', [0.1, 0.86 0.1, 0.1], 'String', 'worst case', 'Color', [0 0 0],"EdgeColor","none","fontsize",9,"Interpreter","latex");
xlim([0 0.8]); ylim([0.8, 1.002]) 
xlabel('$|| \textbf{\emph b}_i - \textbf{\emph b}_{i(1)} ||_2$', 'Interpreter', 'latex');
set(gca,"XTick",0:0.2:0.8); set(gca,"XTickLabel",string(0:0.2:0.8))
set(gca,"YTickLabel",""); set(gca,"YTick",0:0.05:1);
grid on


% Legend on top of the figure
axleg = axes("Position",[0.00 0.82 0.2 0.1]);
plot(0,0,'o',"color", [0.4660 0.6740 0.1880]); hold on
plot(0,1,'sq',"color",[0.4660 0.6740 0.1880]);
ylim([-0.4 1.4])
axis off


% background
axback = axes("Position",[0 0.02 0.5 0.92]);
fill([0 1 1 0],[0 0 1 1],'b',"FaceAlpha",0.08,"EdgeColor","none");
axis off; hold on

axback = axes("Position",[0.5 0.02 0.5 0.92]);
fill([0 1 1 0],[0 0 1 1],'r',"FaceAlpha",0.08,"EdgeColor","none");
axis off; hold on


% Reconstruction errors case 1
i = 1;
ax112 = axes("Position",[0.255 0.58 0.28 0.19]);
imagesc(ax112, X_new(1,:),Y_new(:,1),M(:,:,i)); hold on; caxis([-1.2 1.2]); colormap(ax112,"gray");
fill(xCyl1,yCyl1,'w',xCyl,yCyl2,'w',xCyl,yCyl3,'w','FaceAlpha',1,"EdgeColor","none"); box on
set(gca,"XTickLabel",""); set(gca,"YTickLabel","")
axis equal tight;box on;
xlim([-5 20]); ylim([-5 5]);
title("true")

ax111 = axes("Position",[0.255 0.58 0.28 0.19]);
Int = sqrt(reshape(U_t(1:125561,ind11),size(X_new,1), size(X_new,2)).^2 + reshape(U_t(125562:end,ind11),size(X_new,1), size(X_new,2)).^2);
im = imagesc(ax111, X_new(1,:),Y_new(:,1),Int); hold on; caxis([0 2]); colormap(ax111,"jet"); im.AlphaData = transp; 
fill(xCyl1,yCyl1,'w',xCyl,yCyl2,'w',xCyl,yCyl3,'w','FaceAlpha',1,"EdgeColor","none"); box on
set(gca,"XTickLabel",""); set(gca,"YTickLabel","")
axis equal tight; box on;
xlim([-5 20]); ylim([-5 5]);



i = 2;
ax122 = axes("Position",[0.255 0.29 0.28 0.19]);
im = imagesc(ax122, X_new(1,:),Y_new(:,1),M(:,:,i)); hold on; caxis([-1.2 1.2]); colormap(ax122,"gray");
fill(xCyl1,yCyl1,'w',xCyl,yCyl2,'w',xCyl,yCyl3,'w','FaceAlpha',1,"EdgeColor","none"); box on
set(gca,"XTickLabel",""); set(gca,"YTickLabel","")
axis equal tight;box on;
xlim([-5 20]); ylim([-5 5]);
title("estimation")


ax121 = axes("Position",[0.255 0.29 0.28 0.19]);
Int = sqrt(reshape(U_hat_1(1:125561,ind11),size(X_new,1), size(X_new,2)).^2 + reshape(U_hat_1(125562:end,ind11),size(X_new,1), size(X_new,2)).^2);
im = imagesc(X_new(1,:),Y_new(:,1),Int); hold on; caxis([0 2]); colormap jet; alpha(im,transp)
fill(xCyl1,yCyl1,'w',xCyl,yCyl2,'w',xCyl,yCyl3,'w','FaceAlpha',1,"EdgeColor","none"); box on
set(gca,"XTickLabel","")
set(gca,"YTickLabel","")
axis equal tight
xlim([-5 20]); ylim([-5 5]);
box on;



% Reconstruction errors case 2
i = 5;
ax312 = axes("Position",[0.47 0.58 0.28 0.19]);
im = imagesc(ax312, X_new(1,:),Y_new(:,1),M(:,:,i)); hold on; caxis([-1.2 1.2]); colormap(ax312,"gray");
fill(xCyl1,yCyl1,'w',xCyl,yCyl2,'w',xCyl,yCyl3,'w','FaceAlpha',1,"EdgeColor","none"); box on
set(gca,"XTickLabel","")
set(gca,"YTickLabel","")
axis equal tight
xlim([-5 20]); ylim([-5 5]);
box on;
title("true")

ax311 = axes("Position",[0.47 0.58 0.28 0.19]);
Int = sqrt(reshape(U_t(1:125561,ind12),size(X_new,1), size(X_new,2)).^2 + reshape(U_t(125562:end,ind12),size(X_new,1), size(X_new,2)).^2);
im = imagesc(X_new(1,:),Y_new(:,1),Int); hold on; caxis([0 2]); colormap jet; alpha(im,transp)
fill(xCyl1,yCyl1,'w',xCyl,yCyl2,'w',xCyl,yCyl3,'w','FaceAlpha',1,"EdgeColor","none"); box on
set(gca,"XTickLabel","")
set(gca,"YTickLabel","")
axis equal tight
xlim([-5 20]); ylim([-5 5]);

i = 6;
ax322 = axes("Position",[0.47 0.29 0.28 0.19]);
im = imagesc(ax322, X_new(1,:),Y_new(:,1),M(:,:,i)); hold on; caxis([-1.2 1.2]); colormap(ax322,"gray");
fill(xCyl1,yCyl1,'w',xCyl,yCyl2,'w',xCyl,yCyl3,'w','FaceAlpha',1,"EdgeColor","none"); box on
set(gca,"XTickLabel","")
set(gca,"YTickLabel","")
axis equal tight
xlim([-5 20]); ylim([-5 5]);
box on;
title("estimation")

ax321 = axes("Position",[0.47 0.29 0.28 0.19]);
Int = sqrt(reshape(U_hat_2(1:125561,ind12),size(X_new,1), size(X_new,2)).^2 + reshape(U_hat_2(125562:end,ind12),size(X_new,1), size(X_new,2)).^2);
im = imagesc(X_new(1,:),Y_new(:,1),Int); hold on; caxis([0 2]); colormap jet; alpha(im,transp)
fill(xCyl1,yCyl1,'w',xCyl,yCyl2,'w',xCyl,yCyl3,'w','FaceAlpha',1,"EdgeColor","none"); box on
set(gca,"XTickLabel","")
set(gca,"YTickLabel","")
axis equal tight
xlim([-5 20]); ylim([-5 5]);
box on;



%%%

axcolorbar = axes("Position",[0.35 0.2 1-0.35*2+0.004 0.06]);
colormap("jet"); cbar = colorbar; caxis([0 2])
cbar.Location = "north";
annotation('textbox', [0.46, 0.02 0.1, 0.1], 'String', " $\|\textbf{\emph u}(\textbf{\emph x})\|$", 'Color', [0 0 0],"EdgeColor","none","fontsize",7,"Interpreter","latex");

colormap(ax112,"gray");colormap(ax122,"gray");
colormap(ax312,"gray");colormap(ax322,"gray");

linkaxes([ax111, ax112], 'xy');
linkaxes([ax121, ax122], 'xy');

linkaxes([ax311, ax312], 'xy');
linkaxes([ax321, ax322], 'xy');

ax111.Visible = 'off'; ax111.XTick = []; ax111.YTick = [];
ax121.Visible = 'off'; ax121.XTick = []; ax121.YTick = [];

ax311.Visible = 'off'; ax311.XTick = []; ax311.YTick = [];
ax321.Visible = 'off'; ax321.XTick = []; ax321.YTick = [];
axcolorbar.Visible = 'off'; axcolorbar.XTick = []; axcolorbar.YTick = [];

set(ax111, 'YDir', 'normal');
set(ax121, 'YDir', 'normal');
set(ax112, 'YDir', 'normal');
set(ax122, 'YDir', 'normal');
set(ax311, 'YDir', 'normal');
set(ax321, 'YDir', 'normal');
set(ax312, 'YDir', 'normal');
set(ax322, 'YDir', 'normal');


if saveFig
    if strcmp(DecodingMethod,"MLP"); Name = "Fig6"; else; Name = "Fig7"; end
    exportgraphics(hf , Fig_path + Name +".pdf",'BackgroundColor','none','ContentType',"vector")
end
