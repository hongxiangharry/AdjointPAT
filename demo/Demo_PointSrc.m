% MTR ex1 2016-11-2, code : HXL
% Photoacoustic imaging based on cross-correlation
% Demo for ex1
% 1. Adjoint method
% 2. Modified Time-reversal method
% 3. MTR with low-pass filter
% 4. CINT-like imaging function

clear all; close all;
addpath('../utils');

%% Configuration
saveDataFlag = true; % save all data in the workspace.
Nx_recon = 512; % grid number of ROI per side

%% generate the photoacoustic wavefield

% =========================================================================
% PA FORWARD SIMULATION
% =========================================================================

% load the configuration for the PA forward process
PAconf_PointSrc; 
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% % add noise to the recorded sensor data
signal_to_noise_ratio = inf;
% signal_to_noise_ratio = 40;	% [dB] eqivalent to noise level of 1%
% signal_to_noise_ratio = 20;	% [dB] eqivalent to noise level of 10%
sensor_data = addNoise(sensor_data, signal_to_noise_ratio, 'peak');

%% Adjoint methods

% =========================================================================
% PA ADJOINT RECONSTRUCTION
% =========================================================================

% 1. Time Reversal
Phi_TR = TimeReversal(gather(sensor_data), kgrid.t_array,...
    Nx_recon, x/2, sensor, medium.sound_speed, 'off');

% 2. Back-Projection
Phi_BP =  BackProjection(gather(sensor_data), kgrid.t_array,...
    Nx_recon, x/2, sensor, medium.sound_speed, 'off');

% 3. Truncated Back-Projection
mu = 200/kgrid.t_array(end); % truncate the frequency at the index of 200.
Phi_TBP = TruncatedBackProjection(sensor_data, kgrid.t_array(end),...
    Nx_recon, x/2, sensor, medium.sound_speed, mu, 'on', 'on');

	
%% Normalization and visualization

Ny_recon = Nx_recon;
kgrid_recon = kWaveGrid(Nx_recon, x/Nx_recon, Ny_recon, y/Ny_recon);
disc_x_coord = 0;
disc_y_coord = -R_ROI;
close all;


% TR image
figure1 = figure('position', [500, 500, 1200, 300]);
subaxis(1,4,1, 'MarginTop', 0.02, 'MarginLeft', 0.02, 'MarginRight', 0.02,...
    'Spacing', 0.02, 'MarginBottom',0.02);
hndl = imagesc(kgrid_recon.y_vec*1e3, kgrid_recon.x_vec*1e3, Phi_TR, [-1 1]);
colormap('hot'); axis off; hold on;
scatter_size = 200;
scatter(cart_sensor_mask(2)*1e3, cart_sensor_mask(1)*1e3, scatter_size, 'k^', 'filled');
hold on;
scatter(disc_y_coord*1e3, disc_x_coord*1e3, scatter_size, 'bp', 'filled');
hold on;

% BP image
subaxis(1,4,2, 'MarginTop', 0.02, 'MarginLeft', 0.02, 'MarginRight', 0.02,...
    'Spacing', 0.02, 'MarginBottom',0.02);
hndl = imagesc(kgrid_recon.y_vec*1e3, kgrid_recon.x_vec*1e3, Phi_BP, [-1 1]);
colormap('hot'); axis off; hold on;
scatter_size = 200;
scatter(cart_sensor_mask(2)*1e3, cart_sensor_mask(1)*1e3, scatter_size, 'k^', 'filled');
hold on;
scatter(disc_y_coord*1e3, disc_x_coord*1e3, scatter_size, 'bp', 'filled');
hold on;

% TBP image
subaxis(1,4,3,'MarginTop', 0.02, 'MarginLeft', 0.02, 'MarginRight', 0.02,...
    'Spacing', 0.02, 'MarginBottom',0.02);
hndl = imagesc(kgrid_recon.y_vec*1e3, kgrid_recon.x_vec*1e3, Phi_TBP, [-1 1]);
colormap('hot'); axis off; hold on;
scatter_size = 200;
scatter(cart_sensor_mask(2)*1e3, cart_sensor_mask(1)*1e3, scatter_size, 'k^', 'filled');
hold on;
scatter(disc_y_coord*1e3, disc_x_coord*1e3, scatter_size, 'bp', 'filled'); 
hold on;

% plot a profile for comparison for 3 plot
% [m] location of the slice from the top of reconstruction grid
slice_pos_recon = kgrid_recon.x_size/2;  
slice_pos = kgrid.x_size/2;
subaxis(1,4,4, 'PaddingLeft', 0.05, 'PaddingBottom',0.07);
x_range = round(slice_pos/kgrid.dx/2):round(slice_pos/kgrid.dx);
x_range_recon = round(slice_pos_recon/kgrid_recon.dx/2):round(slice_pos_recon/kgrid_recon.dx);
hndl = plot(kgrid.y_vec(x_range)*1e3, p0(round(slice_pos/kgrid.dx), x_range)/disc_magnitude, 'k--', ...
kgrid_recon.y_vec(x_range_recon)*1e3, Phi_TR(round(slice_pos_recon/kgrid_recon.dx), x_range_recon), 'r-', ...
kgrid_recon.y_vec(x_range_recon)*1e3, Phi_BP(round(slice_pos_recon/kgrid_recon.dx), x_range_recon), 'b-', ...
kgrid_recon.y_vec(x_range_recon)*1e3, Phi_TBP(round(slice_pos_recon/kgrid_recon.dx), x_range_recon), 'g-');
hndl(1).LineWidth = 1.5;
hndl(2).LineWidth = 1.5;
hndl(3).LineWidth = 1.5;
hndl(4).LineWidth = 1.5;
xlabel('y-axis (mm)'); ylabel('Normalized Pressure');
set(gca,'fontsize', 14); leg5 = legend('IP', 'TR', 'BP', 'TBP');
axis tight;
% ax5 = gca;
% outerpos = ax5.OuterPosition;
% ti = ax5.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax5.Position = [left bottom ax_width ax_height];
set(gca, 'YLim', [-1 1.4]);
dataPath = '../output/PointSrc';
if ~exist(dataPath, 'dir')
    mkdir(dataPath);
end
saveas(hndl(1), [dataPath, '/merge_', num2str(num_sensor_points), 'src_snr',...
    num2str(signal_to_noise_ratio), '_N', num2str(Nx), '.png']);
savefig(figure1, [dataPath, '/merge.fig']);
if saveDataFlag == true
    save([dataPath,'/data'], '-v7.3');
end