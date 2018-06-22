% PAconf_OABreast

% Descrption: Configuration of photoacoustic forward process with OA-breast
% phantom model.
% Author: Hongxiang Lin
% Affiliation: Dept. Mechanical Engineering, the University of Tokyo
% 
% Version: 1.0
% Date: 2018-06-20

% load the initial pressure distribution from an image and scale
p0_magnitude = 2;
% p0 = p0_magnitude*loadImage('EXAMPLE_source_two.bmp'); % k-wave example
load model/bloodVesselXY;
p0 = p0_magnitude*phan_proj; % OA-breast model

% assign the grid size and create the computational grid
PML_size = 20;          % size of the PML in grid points
Nx = 1024 - 2*PML_size;  % number of grid points in the x (row) direction
Ny = 1024 - 2*PML_size;  % number of grid points in the y (column) direction
x = 220e-3;              % total grid size [m]
y = 220e-3;              % total grid size [m]
dx = x/Nx;              % grid point spacing in the x direction [m]
dy = y/Ny;              % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% resize the input image to the desired number of grid points
p0_temp = zeros(2*size(p0));
p0_temp(1:size(p0,1),1:size(p0,2)) = p0;
p0 = circshift(p0_temp, floor(size(p0)/2));
p0 = resize(p0, [Nx, Ny]);

% smooth the initial pressure distribution and restore the magnitude
p0 = smooth(kgrid, p0, true);

% define a centered Cartesian circular sensor
sensor_radius = 100e-3;     % [m]
sensor_angle = 2*pi;      % [rad]
sensor_pos = [0, 0];        % [m]
num_sensor_points = 32;
cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points, ...
    sensor_pos, sensor_angle);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s] homogeneous

% assign to sensor structure
sensor.mask = cart_sensor_mask;

% create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% assign to the source structure
source.p0 = p0;

% set the input options
input_args = {'Smooth', false, 'PMLInside', false, 'PlotPML', false, ...
    'PlotSim', true, 'DataCast', 'gpuArray-single'};