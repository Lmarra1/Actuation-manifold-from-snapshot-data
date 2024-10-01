% =============================================================================
% Description: Code of the article "Actuation manifold from snapshots data"
% Authors: Luigi Marra, Guy Y. Cornejo Maceda, Andrea Meilán-Vila, Vanesa Guerrero,
% Salma Rashwan, Bernd R. Noack, Stefano Discetti, and Andrea Ianiro.
% DOI: https://doi.org/10.1017/jfm.2024.593
% Dataset DOI: 10.5281/zenodo.12802191.
% GitHub: https://github.com/Lmarra1/Actuation-manifold-from-snapshot-data.git
% Description:
% This script performs the plot of the results of ISOMAP dimensionality reduction on a
% dataset of snapshots of the fluidic pinball with control. Plots of the
% encoding step.
% =============================================================================

%% Clear previous data and settings
close all; clear; clc; warning off; % Disable warnings

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

%% Load data for plotting
% Update this path to the location where the data downloaded from Zenodo (DOI: [To be received]) is stored
Data_path = "SPECIFY YOUR DATA PATH";

Fig_path = "ArticleFigures";
if ~exist(Fig_path, 'dir')
    mkdir(Fig_path)
end

if ~exist("toolbox_image",'dir')
    disp("To plot the LIC of the pseudomodes you need to download the toolbox_image:")
    disp("Gabriel Peyre (2024). Toolbox image" +...
        "(https://www.mathworks.com/matlabcentral/fileexchange/16201-toolbox-image),"+...
        "MATLAB Central File Exchange.")
else
    addpath(genpath("toolbox_image"));
end

ReadH5(fullfile(Data_path, "TrainingDataset_1.h5"));          % Load the training dataset data
ReadH5(fullfile(Data_path, "IsomapVarke.h5"));                % Load the initial data needed for optimal k_e selection
ReadH5(fullfile(Data_path, "IsomapResults.h5"));              % Load the Isomap results for k_e = 40
ReadH5(fullfile(Data_path, "Grid.h5"));                       % Load the grid data


%% Other general variables
% Define circumference and coordinates for plotting the fluidic pinball cylinders
theta = 2*pi*(0:100)/100;           % Define theta for circumference
xCyl1 = -1.299 + 0.5*cos(theta);    % X coordinate for cylinder 1
yCyl1 = 0.5*sin(theta);             % Y coordinate for cylinder 1
xCyl  = 0.5*cos(theta);             % X coordinate for cylinder 2 and 3
yCyl2 = 0.75 + 0.5*sin(theta);      % Y coordinate for cylinder 2
yCyl3 = -0.75 + 0.5*sin(theta);     % Y coordinate for cylinder 3

% Define limits for plotting low-dim embedding
glim = [-24 12; repmat([-12 12],9,1)];

% Extract kiki parameters for plotting
BT      = p(1,:); % Data for boat tailinclc
% g
MAGNUS  = p(2,:); % Data for MAGNUS
SP      = p(3,:); % Data for forward stagnation point



%% Plots of figure n*2: Manifold sections and k_e selection
saveFig = false; % Flag to save the figure

% Fig2a: Frobenius norm of the geodesic distance matrix as k_e varies
figure('units','centimeters','position',[10, 3, width, width*9/(16)*1.45]); hold on
plot(ke, FrobNormDg,"k.-","MarkerSize",5,'linewidth',0.5);
xl = xline(40,"k--","linewidth",0.5); xl.Color = [0 0 0 0.5];
xlabel("$k_e$"); ylabel("$||\mathbf{D}_G||_F$"); grid on
fill([0 40 40 0],[-1000000 -1000000 1000000 1000000],"b","FaceAlpha",0.07,"EdgeColor","none")
xlim([0 150]); ylim([70000 100000]);
set(gca,"fontsize",fs); set(gca,'ytick',70000:10000:100000);
plot(ke(1), FrobNormDg(1),"r.","MarkerSize",8);

if saveFig
    Name = "Fig2a";
    exportgraphics(hf , Fig_path + Name + ".pdf",'BackgroundColor','none','ContentType',"vector")
end


% Fig2b: Number of points connected in the KNN-graph
hf = figure('units','centimeters','position',[10, 3, width, width*9/(16)*1.45]); hold on
plot(1:60, NpointsConn(1:60)/6860*100,"k.-","MarkerSize",5,'linewidth',0.5);
xl = xline(40,"k--","linewidth",0.5); xl.Color = [0 0 0 0.5];
xlabel("$k_e$"); ylabel("$\%$"); grid on
fill([0 40 40 0],[-1000000 -1000000 1000000 1000000],"b","FaceAlpha",0.07,"EdgeColor","none")
xlim([0 60]); ylim([0, 100])
set(gca,"fontsize",fs); set(gca,'ytick',0:25:100); set(gca,'xtick',0:20:60);

if saveFig
    Name = "Fig2b";
    exportgraphics(hf , Fig_path + Name + ".pdf",'BackgroundColor','none','ContentType',"vector")
end


% Fig2c: Residual variance
hf = figure('units','centimeters','position',[10, 3, width, width*9/(16)*1.45]); hold on
plot(1:10, Rv,"k.-","MarkerSize",20,'linewidth',1);
xlabel("$n$"); ylabel("$R_v$"); grid on
xlim([1 10]); ylim([0, 1])
set(gca,"fontsize",fs); set(gca,'xtick',1:10); set(gca,'ytick',0:0.2:1);

if saveFig
    Name = "Fig2c";
    exportgraphics(hf , Fig_path + Name + ".pdf",'BackgroundColor','none','ContentType',"vector")
end


% Figd: All manifold sections
Gamma2plot_vec = [1 2; 1 3; 2 3; 1 4; 2 4; 3 4; 1 5; 2 5; 3 5; 4 5];

for i = 1:size(Gamma2plot_vec,1)
    figure('units','centimeters','position',[10, 3, width, width*9/(16)]);
    Gamma2plot   = Gamma2plot_vec(i,:);
    scatter(Gamma(Gamma2plot(1),:), Gamma(Gamma2plot(2),:),PointSize ,"blue", "Marker",".");
    box on; grid on; hold on; axis equal
    xlabel("$\gamma_" + string(Gamma2plot_vec(i,1)) + "$")
    ylabel("$\gamma_" + string(Gamma2plot_vec(i,2)) + "$")
    xlim(glim(Gamma2plot(1),:)); ylim(glim(Gamma2plot(2),:));

    if saveFig
        Name = "Figd_" + string(i);
        exportgraphics(hf , Fig_path + Name + ".pdf",'BackgroundColor','none','ContentType',"vector")
    end

end



%% Plots of figure n*3: Horizontal Symmetry Interpretation
saveFig = false; % Flag to save the figure


% Fig3b: plot of the manifold section gamma_1-gamma_2
Gamma2plot   = [1,2];
saveFig = false;
figure('units','centimeters','position',[10, 3, width, width*9/(16)*1.2]);
sc = scatter(Gamma(Gamma2plot(1), :), Gamma(Gamma2plot(2), :), 8, [0 0 0], 'Marker', 'o');
axis equal; axis off; axis tight;
xlim(glim(Gamma2plot(1),:)); ylim(glim(Gamma2plot(2),:));
sc.MarkerFaceColor = [0 0 0];
sc.MarkerEdgeColor = "none";

if saveFig
    Name = "Fig3b";
    exportgraphics(hf , Fig_path + Name + ".pdf",'BackgroundColor','none','ContentType',"vector")
end



% Fig3a: plot of the manifold section gamma_1-gamma_2 color coded with b_1, b_2, b_3
figure('units','centimeters','position',[10, 3, width, width*9/(16*4)*2.8]);
Gamma2plot   = [1,2];


PLOT = tiledlayout(1, 3);

nexttile
scatter(Gamma(Gamma2plot(1),:), Gamma(Gamma2plot(2),:), PointSize, b(1,:) ,"Marker",".");
xlm = glim(Gamma2plot(1),:); ylm = glim(Gamma2plot(2),:);
xlim(xlm); ylim(ylm)
set(gca,'xtick',xlm(1):6:xlm(2)); set(gca,'xticklabel',xlm(1):6:xlm(2));
set(gca,'ytick',ylm(1):6:ylm(2)); set(gca,'yticklabel',ylm(1):6:ylm(2));
ylabel("$\gamma_"+ string(Gamma2plot(2)) + "$");
axis equal; xlim(xlm); ylim(ylm);
colormap(jet(7))
caxis([-3.5 3.5])
title("$b_1$","FontSize",fs)


nexttile
scatter(Gamma(Gamma2plot(1),:), Gamma(Gamma2plot(2),:), PointSize, b(2,:) ,"Marker",".");
xlabel("$\gamma_"+ string(Gamma2plot(1)) + "$");
axis equal; xlim(xlm); ylim(ylm)
colormap(jet(7))
caxis([-3.5 3.5])
title("$b_2$","FontSize",fs)
set(gca,'xtick',xlm(1):6:xlm(2)); set(gca,'xticklabel',xlm(1):6:xlm(2));
set(gca,'ytick',ylm(1):6:ylm(2));
set(gca,'yticklabel',"");


nexttile
scatter(Gamma(Gamma2plot(1),:), Gamma(Gamma2plot(2),:), PointSize, b(3,:) ,"Marker",".");
axis equal; xlim(xlm); ylim(ylm)
colormap(jet(7))
caxis([-3.5 3.5])
title("$b_3$","FontSize",fs)
set(gca,'xtick',xlm(1):6:xlm(2)); set(gca,'xticklabel',xlm(1):6:xlm(2));
set(gca,'ytick',ylm(1):6:ylm(2));
set(gca,'yticklabel',"");
cbh = colorbar;
cbh.Ticks = -3:3;
cbh.TickLabels = string(-3:3);
cbh.Layout.Tile = 'north';
cbh.Position = cbh.Position + [0 0 0 -0.02];



PLOT.TileSpacing = 'compact';
PLOT.Padding     = 'compact';

if saveFig
    Name = "Fig3a";
    exportgraphics(hf , Fig_path + Name + ".pdf",'BackgroundColor','none','ContentType',"vector")
end





%% Plots of figure n*4: Interpretation of manifold coordinates
saveFig = false;

hf = figure('units','centimeters','position',[10, 3, 13.2, 7.2]);
PLOT = tiledlayout(2,3);

% 1 1
nexttile
Gamma2plot   = [1,3,4];
contour2show = "BT";

scatter3(Gamma(Gamma2plot(1),:), Gamma(Gamma2plot(2),:), Gamma(Gamma2plot(3),:), PointSize, eval(contour2show),"Marker",".");
xlm = glim(Gamma2plot(1),:); ylm = glim(Gamma2plot(2),:); zlm = glim(Gamma2plot(3),:);
axis equal; xlim(xlm); ylim(ylm); zlim(zlm);
set(gca,'xtick',xlm(1):6:xlm(2)); set(gca,'xticklabel',xlm(1):6:xlm(2));
set(gca,'ytick',ylm(1):6:ylm(2)); set(gca,'yticklabel',ylm(1):6:ylm(2));
set(gca,'ztick',zlm(1):6:zlm(2)); set(gca,'zticklabel',zlm(1):6:zlm(2));
xlabel("$\gamma_"+ string(Gamma2plot(1)) + "$");
yl = ylabel("$\gamma_"+ string(Gamma2plot(2)) + "$");
zlabel("$\gamma_"+ string(Gamma2plot(3)) + "$");
view([-15 20]); colormap(jet(16))
yl_pos = get(yl, 'Position');
yl_pos = yl_pos + [5, 15, 3]; % Adjust the values as needed
set(yl, 'Position', yl_pos);
caxis([-3 3])
cbh(1) = colorbar;
cbh(1).Ticks = -3:3;
cbh(1).TickLabels = string(-3:3);
cbh(1).Location = "eastoutside";
cbh(1).Position = cbh(1).Position + [-0.05 0.0145 0 0];
title(cbh(1),"$p_1$","interpreter","latex","FontSize",fs)
cbh(1).Position = (cbh(1).Position + [0.025 0.0 0 0]).*[1 1 0.4 1];
cbh(1).Position(4) = 0.25;
cbh(1).Position(2) = 0.68;


% 1 2
nexttile
contour2show = "MAGNUS";
Gamma2plot   = [2,3,4];
scatter3(Gamma(Gamma2plot(1),:), Gamma(Gamma2plot(2),:), Gamma(Gamma2plot(3),:), PointSize, eval(contour2show),"Marker",".");
xlm = glim(Gamma2plot(1),:); ylm = glim(Gamma2plot(2),:); zlm = glim(Gamma2plot(3),:);
axis equal;xlim(xlm); ylim(ylm); zlim(zlm);
set(gca,'xtick',xlm(1):6:xlm(2)); set(gca,'xticklabel',xlm(1):6:xlm(2));
set(gca,'ytick',ylm(1):6:ylm(2)); set(gca,'yticklabel',ylm(1):6:ylm(2));
set(gca,'ztick',zlm(1):6:zlm(2)); set(gca,'zticklabel',zlm(1):6:zlm(2));
xl = xlabel("$\gamma_"+ string(Gamma2plot(1)) + "$");
yl = ylabel("$\gamma_"+ string(Gamma2plot(2)) + "$");
zlabel("$\gamma_"+ string(Gamma2plot(3)) + "$");
view([-15 20]);
colormap(jet(16))
Labelshift = [2,  1,  7];
yl_pos = get(yl, 'Position');
yl_pos = yl_pos + Labelshift; % Adjust the values as needed
set(yl, 'Position', yl_pos);
Labelshift = [0,  4,  1];
xl_pos = get(xl, 'Position');
xl_pos = xl_pos + Labelshift; % Adjust the values as needed
set(xl, 'Position', xl_pos);
caxis([-9 9])
cbh(1) = colorbar;
cbh(1).Ticks = -9:3:9;
cbh(1).TickLabels = string(-9:3:9);
title(cbh(1),"$p_2$","interpreter","latex","FontSize",fs)
cbh(1).Location = "Eastoutside";
cbh(1).Position = (cbh(1).Position + [-0.01 0.0 0 0]).*[1 1 0.4 1];
cbh(1).Position(4) = 0.25;
cbh(1).Position(2) = 0.68;


% 1 3
nexttile
Gamma2plot   = [5,3,4];
contour2show = "MAGNUS";
scatter3(Gamma(Gamma2plot(1),:), Gamma(Gamma2plot(2),:), Gamma(Gamma2plot(3),:), PointSize, eval(contour2show),"Marker",".");
xlm = glim(Gamma2plot(1),:); ylm = glim(Gamma2plot(2),:); zlm = glim(Gamma2plot(3),:);
axis equal;xlim(xlm); ylim(ylm); zlim(zlm);
set(gca,'xtick',xlm(1):6:xlm(2)); set(gca,'xticklabel',xlm(1):6:xlm(2));
set(gca,'ytick',ylm(1):6:ylm(2)); set(gca,'yticklabel',ylm(1):6:ylm(2));
set(gca,'ztick',zlm(1):6:zlm(2)); set(gca,'zticklabel',zlm(1):6:zlm(2));
xl = xlabel("$\gamma_"+ string(Gamma2plot(1)) + "$");
yl = ylabel("$\gamma_"+ string(Gamma2plot(2)) + "$");
zlabel("$\gamma_"+ string(Gamma2plot(3)) + "$");
view([-15 20]); colormap(jet(16))
yl_pos = get(yl, 'Position');
yl_pos = yl_pos + [2, 0, 7];
set(yl, 'Position', yl_pos);
xl_pos = get(xl, 'Position');
xl_pos = xl_pos + [0,  3,  2];
set(xl, 'Position', xl_pos);
caxis([-3 3])
cbh(1) = colorbar;
cbh(1).Ticks = -3:3;
cbh(1).TickLabels = string(-3:3);
cbh(1).Location = "Eastoutside";
cbh(1).Position = (cbh(1).Position + [+0.04 0.0 0 0]).*[1 1 0.4 1];
cbh(1).Position(4) = 0.25;
cbh(1).Position(2) = 0.68;
title(cbh(1),"$p_3$","interpreter","latex","FontSize",fs)



% 2 1
pts = 1;
nexttile
scatter(BT,Cd,pts,"Marker",".");
xlabel("$p_1$"); ylabel("$C_D$")
ylim([0 13]); xlim([-4 4])
set(gca,'xtick',-4:2:4); set(gca,'xticklabel',string(-4:2:4));
set(gca,'ytick',0:4:12); set(gca,'yticklabel',string(0:4:12));
axis square;grid on


% 2 2
nexttile
scatter(MAGNUS,Cl,pts,"Marker",".");

xlabel("$p_2$"); ylabel("$C_L$")
set(gca,'xtick',-10:5:10); set(gca,'xticklabel',string(-10:5:10));
set(gca,'ytick',-20:10:20); set(gca,'yticklabel',string(-20:10:20));
xlim([-10 10]); ylim([-20 20]);
axis square;grid on


% 2 3
nexttile
Gamma2plot   = [5,3,4];
contour2show = "Cl";
scatter3(Gamma(Gamma2plot(1),:), Gamma(Gamma2plot(2),:), Gamma(Gamma2plot(3),:), PointSize, eval(contour2show),"Marker",".");
xlm = glim(Gamma2plot(1),:); ylm = glim(Gamma2plot(2),:); zlm = glim(Gamma2plot(3),:);
axis equal;xlim(xlm); ylim(ylm); zlim(zlm);
set(gca,'xtick',xlm(1):6:xlm(2)); set(gca,'xticklabel',xlm(1):6:xlm(2));
set(gca,'ytick',ylm(1):6:ylm(2)); set(gca,'yticklabel',ylm(1):6:ylm(2));
set(gca,'ztick',zlm(1):6:zlm(2)); set(gca,'zticklabel',zlm(1):6:zlm(2));
xlabel("$\gamma_"+ string(Gamma2plot(1)) + "$");
yl = ylabel("$\gamma_"+ string(Gamma2plot(2)) + "$");
zlabel("$\gamma_"+ string(Gamma2plot(3)) + "$");
view([-15 20]); colormap(jet(16))
yl_pos = get(yl, 'Position');
yl_pos = yl_pos + [3, 3, 6]; % Adjust the values as needed
set(yl, 'Position', yl_pos);
caxis([-18 18])
cbh(1) = colorbar;
cbh(1).Ticks = -18:9:18;
cbh(1).TickLabels = string(-18:9:18);
cbh(1).Location = "Eastoutside";
cbh(1).Position = (cbh(1).Position + [+0.05 0.0 0 0]).*[1 1 0.4 1];
cbh(1).Position(4) = 0.25;
cbh(1).Position(2) = 0.18;
% annotation('textbox', [0.822 0.04 .001 .01], 'String', '$C_L$', 'Interpreter', 'latex', 'FontSize', fs-1, 'EdgeColor', 'none','FitBoxToText','on');
title(cbh(1),"$C_L$","interpreter","latex","FontSize",fs)

PLOT.TileSpacing = 'compact';
PLOT.Padding     = 'compact';

if saveFig
    Name = "Fig4";
    exportgraphics(hf , Fig_path + Name + ".pdf",'BackgroundColor','none','ContentType',"vector")
end







%% Plots of figure n*5: Pseudomodes LIC

% Perform pseudomodes LIC
saveFig = false;

cmap = jet(16);
xsq = linspace(-5,20,300);
ysq = linspace(-5,5 ,300);
[Xsq,Ysq] = meshgrid(xsq,ysq);
Phi_u(logical(Mask(:)),:) = 1;
Phi_v(logical(Mask(:)),:) = 1;

for i = 1:5

    Ut = reshape(Phi_u(:,i),size(X_new,1),size(X_new,2));
    Vt = reshape(Phi_v(:,i),size(X_new,1),size(X_new,2));

    Ut = interp2(X_new,Y_new,Ut,Xsq,Ysq,"linear");
    Vt = interp2(X_new,Y_new,Vt,Xsq,Ysq,"linear");

    VelF = cat(3, Vt, Ut);

    VelF = VelF./sqrt(sum(VelF.^2,3));

    % Length of the convolution
    w = 80;

    % Set options
    options = struct();
    options.spot_size = 5;    % Set the size of the features

    % Perform Line Integral Convolution
    M(:,:,i) = perform_lic(VelF, w, options);
    disp("Pseudomode " + string(i) + ": LIC performed")

end

%% Plot the pseudomodes LIC

for i = 1:5

    hf = figure('units','centimeters','position',[10, 1.5, 4, 2]);

    Umi = reshape(Phi_u(:,i),size(X_new,1),size(X_new,2));
    Vmi = reshape(Phi_v(:,i),size(X_new,1),size(X_new,2));
    ax11 = axes("Position",[0.05 0.05 0.9 0.9]);
    pc=imagesc(ax11,Xsq(1,:),Ysq(:,1).',M(:,:,i));
    hold on; shading interp; axis equal tight;box on;
    fill(xCyl1,yCyl1,[1 1 1 ],xCyl,yCyl2,[1 1 1],xCyl,yCyl3,[1 1 1],"EdgeColor","none");
    set(gca,'fontsize',fs);
    set(gca,"ytick",-5:5:5); set(gca,"yticklabel",string(-5:5:5))
    set(gca,"xtick",-5:5:20); set(gca,"xticklabel",string(-5:5:20))
    xlim([-5 20]); ylim([-5 5]);
    colormap(ax11,"gray");
    caxis([-1.2 1.2])
    axis off

    ax12 = axes("Position",[0.05 0.05 0.9 0.9]);
    pc=imagesc(ax12,X_new(1,:),Y_new(:,1).',sqrt(Umi.^2+Vmi.^2));
    hold on; shading interp; axis equal tight;box on;
    fill(xCyl1,yCyl1,[1 1 1],xCyl,yCyl2,[1 1 1],xCyl,yCyl3,[1 1 1],"EdgeColor","none");
    set(gca,'fontsize',fs);
    set(gca,"ytick",-5:5:5); set(gca,"yticklabel",string(-5:5:5))
    set(gca,"xtick",-5:5:20); set(gca,"xticklabel",string(-5:5:20))
    xlim([-5 20]); ylim([-5 5]);
    colormap(ax12,"jet(16)"); caxis([0 20])
    pc.AlphaData = 0.5;
    axis off;

    % Link axes
    linkaxes([ax11, ax12], 'xy');
    ax12.Visible = 'off'; ax12.XTick = []; ax12.YTick = [];
    if saveFig
        Name = "Fig5_" + string(i);
        exportgraphics(hf , Fig_path + Name +".pdf",'BackgroundColor','none','ContentType',"vector")
    end

end



