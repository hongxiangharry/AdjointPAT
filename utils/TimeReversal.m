function [Phi] = TimeReversal(P, t_array, Nd, Rd, sensor, sound_speed, GPUmode)
% TimeReversal
% Description: calculate the initial pressure distribution for PAT by using
% the Time-Reversal method. ROI is a square centered at the origin. The
% code is originated from the method of 2D time reversal reconstruction in
% the k-Wave toolbox.
% Input:
%   P: received waveform data
%   t_array: time sequence for the received waveform data
%   Nd: grid points per side of ROI
%   Rd: half length of a side of ROI
%   sensor: sensor information defined by the k-Wave toolbox
%   sound_speed: sound speed of the background
%   GPUmode: on/off of GPU mode
% Output:
%   Phi: the TR imaging function
%
% Author: Hongxiang Lin
% Affiliation: Dept. Mechanical Engineering, the University of Tokyo
% 
% Version: 1.1
% Date: 2018-06-20
% 
% Original description from the source code
% k-Wave\examples\example_pr_2D_TR_circular_sensor.m
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

if nargin < 7
    GPUmode = 'off';
end
Nx = Nd;               % number of grid points in the x (row) direction
Ny = Nd;               % number of grid points in the y (column) direction
dx = 2*Rd/Nx;          % grid point spacing in the x direction [m]
dy = 2*Rd/Ny;          % grid point spacing in the y direction [m]
kgrid_recon = kWaveGrid(Nx, dx, Ny, dy);
kgrid_recon.t_array = t_array;
medium.sound_speed = sound_speed;
source_recon.p0 = 0; % reset the initial pressure
sensor.time_reversal_boundary_data = P; % assign the time reversal data
if strcmp(GPUmode, 'off') % CPU mode
    input_args = {'Smooth', false, 'PMLInside', false, 'PlotPML', false,...
        'PlotSim', false, 'DisplayMask', cart2grid(kgrid_recon, sensor.mask)};
    % run the time-reversal reconstruction
    Phi = kspaceFirstOrder2D(kgrid_recon, medium, source_recon, sensor, input_args{:});
elseif strcmp(GPUmode, 'on') % GPU mode
    input_args = {'Smooth', false, 'PMLInside', false, 'PlotPML', false,...
        'PlotSim', false, 'DisplayMask', cart2grid(kgrid_recon, sensor.mask),...
        'DataCast', 'gpuArray-single'};
    % run the time-reversal reconstruction
    Phi = kspaceFirstOrder2D(kgrid_recon, medium, source_recon, sensor, input_args{:});
    Phi = gather(Phi);
end
Phi(isnan(Phi)) = 0; % remove the exceptional nan-value.
Phi = real(Phi/max(max(abs(Phi))));