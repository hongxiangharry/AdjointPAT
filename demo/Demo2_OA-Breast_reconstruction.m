% Demo2_OA-Breast_reconstruction
% Description: A demo for the adjoint-method reconstruction of the OA-breast
% phantom. Three adjoint methods for PAT are compared:
%   1. Time-Reversal (TR)
%   2. Back-Projection (BP)
%   3. Truncated Back-Projection (TBP)
% 
% Forward process of synthesizing photoacoustic wave signals are simulated 
% by the k-Wave toolbox. This is downloaded from the website:
%       http://www.k-wave.org/
% Most valuables emerging in this script are directly used as their
% original definitions specified in the k-Wave toolbox, such as sensor,
% source, medium etc. See the details in the link for help:
% http://www.k-wave.org/documentation/kspaceFirstOrder2D.php
% 
% The OA-breast phantom is an open database of breast phantom released by 
% Computational Imaging Science Laboratory, Washington Unversity in
% St.Louis. See the link below for details:
% https://anastasiolab.wustl.edu/downloadable-content/oa-breast-database/
%
% Author: Hongxiang Lin
% Affiliation: Dept. Mechanical Engineering, the University of Tokyo
% 
% Version: 1.0
% Date: 2018-06-21

clear all; close all;
addpath('../utils');

%% Configuration
saveDataFlag = true; % save all data in the workspace.
Nx_recon = 256; % grid number of ROI per side
slice_pos = 100e-3;  % [m] location of the slice from top of grid [m]

%% Generate photoacoustic wavefield using the k-Wave toolbox

% =========================================================================
% PA FORWARD SIMULATION
% =========================================================================

% load the configuration for the PA forward process
PAconf_OABreast; 

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% add noise to the recorded sensor data
signal_to_noise_ratio = inf;
% signal_to_noise_ratio = 40;	% [dB] eqivalent to noise level of 1%
% signal_to_noise_ratio = 20;	% [dB] eqivalent to noise level of 10%
sensor_data = addNoise(sensor_data, signal_to_noise_ratio, 'peak');

%% Adjoint methods

% =========================================================================
% PA ADJOINT RECONSTRUCTION
% =========================================================================

% 1. Time Reversal
Phi_TR = TimeReversal(sensor_data, kgrid.t_array,...
    Nx_recon, x/2, sensor, medium.sound_speed, 'on');

% 2. Back-Projection
Phi_BP =  BackProjection(sensor_data, kgrid.t_array,...
    Nx_recon, x/2, sensor, medium.sound_speed, 'on');

% 3. Truncated Back-Projection
mu = 200/kgrid.t_array(end); % truncate the frequency at the index of 200.
Phi_TBP = TruncatedBackProjection(sensor_data, kgrid.t_array(end),...
    Nx_recon, x/2, sensor, medium.sound_speed, mu, 'on', 'on');

%% Normalization and visualization

Ny_recon = Nx_recon;
x_ROI = x/2;
y_ROI = y/2;
x_idx = floor((x-x_ROI)/x*Nx_recon/2);
y_idx = floor((y-y_ROI)/y*Ny_recon/2);
kgrid_recon = kWaveGrid(Nx_recon, x/Nx_recon, Ny_recon, y/Ny_recon);
close all;

% TR image
figure1 = figure('position', [500, 500, 1200, 300]);
subaxis(1,4,1, 'MarginTop', 0.02, 'MarginLeft', 0.02, 'MarginRight', 0.02,...
    'Spacing', 0.02, 'MarginBottom',0.02);
hndl = imagesc(kgrid_recon.y_vec(y_idx:Ny_recon-y_idx)*1e3,...
    kgrid_recon.x_vec(x_idx:Nx_recon-x_idx)*1e3,...
    Phi_TR(x_idx:Nx_recon-x_idx, y_idx:Ny_recon-y_idx), [-1 1]);
colormap('hot'); ylabel('x-axis [mm]'); xlabel('y-axis [mm]');
set(gca,'fontsize', 14);
axis off;

% BP image
subaxis(1,4,2, 'MarginTop', 0.02, 'MarginLeft', 0.02, 'MarginRight', 0.02,...
    'Spacing', 0.02, 'MarginBottom',0.02);
hndl = imagesc(kgrid_recon.y_vec(y_idx:Ny_recon-y_idx)*1e3,...
    kgrid_recon.x_vec(x_idx:Nx_recon-x_idx)*1e3,...
    Phi_BP(x_idx:Nx_recon-x_idx, y_idx:Ny_recon-y_idx), [-1 1]);
colormap('hot'); ylabel('x-axis [mm]'); xlabel('y-axis [mm]');
set(gca,'fontsize', 14);
axis off;

% TBP image
subaxis(1,4,3, 'MarginTop', 0.02, 'MarginLeft', 0.02, 'MarginRight', 0.02,...
    'Spacing', 0.02, 'MarginBottom',0.02);
hndl = imagesc(kgrid_recon.y_vec(y_idx:Ny_recon-y_idx)*1e3,...
    kgrid_recon.x_vec(x_idx:Nx_recon-x_idx)*1e3,...
    Phi_TBP(x_idx:Nx_recon-x_idx, y_idx:Ny_recon-y_idx), [-1 1]);
colormap('hot'); ylabel('x-axis [mm]'); xlabel('y-axis [mm]');
set(gca,'fontsize', 14);
axis off;

% plot a profile for comparison of TR, BP, and TBP
grid_ext = kgrid.y_vec*1e3;
grid_recon = kgrid_recon.y_vec*1e3;
p0_pf_ext = p0(round(slice_pos/kgrid.dx), :)/p0_magnitude;
p0_pf_TR  = Phi_TR(round(slice_pos/kgrid_recon.dx), :);
p0_pf_BP  = Phi_BP(round(slice_pos/kgrid_recon.dx), :);
p0_pf_TBP = Phi_TBP(round(slice_pos/kgrid_recon.dx), :);

subaxis(1,4,4, 'PaddingLeft', 0.05, 'PaddingBottom', 0.07);
hndl = plot(grid_ext, p0_pf_ext, 'k--', ...
            grid_recon, p0_pf_TR, 'r-', ...
            grid_recon, p0_pf_BP, 'b-',...
            grid_recon, p0_pf_TBP, 'g-');
hndl(1).LineWidth = 1.5;
hndl(2).LineWidth = 1.5;
hndl(3).LineWidth = 1.5;
hndl(4).LineWidth = 1.5;
xlabel('y-axis (mm)'); ylabel('Normalized Pressure'); 
set(gca,'fontsize', 14); legend('IP', 'TR', 'BP', 'TBP'); axis tight;
set(gca, 'XLim', [0 75]); set(gca, 'YLim', [-.4 1.2]);

dataPath = '../output/OA-Breast';
if ~exist(dataPath, 'dir')
    mkdir(dataPath);
end
saveas(hndl(1), [dataPath,'/merge_snr', num2str(signal_to_noise_ratio), '.png']);

savefig(figure1, [dataPath, '/merge.fig']);
if saveDataFlag == true
    save([dataPath,'/data'], '-v7.3');
end
