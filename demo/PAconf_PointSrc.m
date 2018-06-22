% PAconf_PointSrc

% Descrption: Configuration of photoacoustic forward process with a point
% source.
% 
% Author: Hongxiang Lin
% Affiliation: Dept. Mechanical Engineering, the University of Tokyo
% 
% Version: 1.0
% Date: 2018-06-20

% assign the grid size and create the computational grid
PML_size = 20;          % size of the PML in grid points
Nx = 1024 - 2*PML_size;  % number of grid points in the x (row) direction
Ny = 1024 - 2*PML_size;  % number of grid points in the y (column) direction
x = 220e-3;              % total grid size [m]
y = 220e-3;              % total grid size [m]
dx = x/Nx;              % grid point spacing in the x direction [m]
dy = y/Ny;              % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define a centered Cartesian circular sensor
sensor_radius = 100e-3;     % [m]
sensor_angle = 2*pi;      % [rad]
sensor_pos = [0, 0];        % [m]
num_sensor_points = 1;
cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points, ...
    sensor_pos, sensor_angle);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s] homogeneous sound speed

% assign to sensor structure
sensor.mask = cart_sensor_mask;

% create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% assign to the source structure

%%%%%%%%%%%%%%%%% single point source %%%%%%%%%%%%%%%%%%%%
disc_magnitude = 50; % [au]
disc_src_num_per_row = 4; % source numbers per row
R_ROI = x/8;
disc_x_array = floor(Nx/x*x/2); % [grid points]
disc_y_array = floor(Ny/y*(y/2-R_ROI)); % [grid points]
disc_radius = floor(Nx/4/disc_src_num_per_row/4);   % [grid points], 1/4 
p0 = zeros(Nx, Ny);
disc = disc_magnitude*makeDisc(Nx, Ny, disc_x_array(1), disc_y_array(1), disc_radius);
p0 = p0+disc;
source.p0 = p0;

% % smooth the initial pressure distribution and restore the magnitude
% FigHandle = figure; 
% % set(FigHandle, 'Position', [100, 100, 1049, 600]);
% % subplot(1,2,1); 
% imagesc(source.p_mask, [-1 1]); axis square; colorbar; colormap('hot');

% set the input options
input_args = {'Smooth', false, 'PMLInside', false, 'PlotPML', false, ...
    'PlotSim', true, 'DataCast', 'gpuArray-single'};