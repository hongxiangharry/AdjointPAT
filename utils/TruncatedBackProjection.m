function [Phi] = TruncatedBackProjection(P, T, Nd, Rd, sensor, sound_speed,...
    mu, GPUmode, PlotMode)
% TruncatedBackProjection
% Description: calculate the initial pressure distribution for PAT by using
% Truncated Back-Projection. ROI is a square centered at the origin.
% Input:
%   P: received waveform data
%   T: time duration
%   Nd: grid points per side of ROI
%   Rd: half length of a side of ROI
%   sensor: sensor information defined by the k-Wave toolbox
%   sound_speed: sound speed of the background
%   mu: highest frequency of truncated band of signals
%   GPUmode: on/off of GPU mode
%   PlotMode: on/off for plotting the summation of imaging function for
%   the frequency components
% Output:
%   Phi: the BP imaging function
% 
% Author: Hongxiang Lin
% Affiliation: Dept. Mechanical Engineering, the University of Tokyo
% 
% Version: 1.0
% Date: 2018-06-20

if nargin < 8
    GPUmode = 'off';
elseif nargin < 9
	PlotMode = 'off';
end
Nx = Nd;           % number of grid points in the x (row) direction
Ny = Nd;           % number of grid points in the y (column) direction
dx = 2*Rd/Nx;          % grid point spacing in the x direction [m]
dy = 2*Rd/Ny;          % grid point spacing in the y direction [m]
kgrid_recon = kWaveGrid(Nx, dx, Ny, dy);
num_sensor_points = size(sensor.mask, 2);
P = fft(P, size(P, 2), 2);
% GPU version:
if strcmp(GPUmode, 'on')
	dGrid = gpuArray.linspace(-Rd, Rd, Nd);
	Phi = zeros(Nd, Nd, 'single', 'gpuArray');
else
	dGrid = linspace(-Rd, Rd, Nd);
	Phi = zeros(Nd, Nd);	
end
fq_no = 0;
if strcmp(PlotMode, 'on')
    figure;
    imagesc(kgrid_recon.y_vec*1e3, kgrid_recon.x_vec*1e3,log(abs((Phi))));
    chndl = colorbar; title(['TBP, index of Fq = ', num2str(fq_no)]);
    colormap('hot'); ylabel('x-axis [mm]'); xlabel('y-axis [mm]');
    title(chndl, 'log'); set(gca,'fontsize', 16); axis square;
    axis tight manual; ax = gca; ax.NextPlot = 'replaceChildren';
    drawnow
    F(1) = getframe;
end

loops = ceil(mu*T+1); % Truncated index of frequency components
if loops < 1 || loops > size(P, 2)
	error('Truncated index should be a positive integer small than number of time steps');
end
F(loops) = struct('cdata',[],'colormap',[]);
for fq_no = 2:loops
    fq = (fq_no-1)/T;
	if strcmp(GPUmode, 'on')
		GreenFun = zeros(Nd, Nd, num_sensor_points, 'single', 'gpuArray');
	else
		GreenFun = zeros(Nd, Nd, num_sensor_points);
	end		
    for idx_sensor = 1:num_sensor_points
        r_sensor = sensor.mask(:, idx_sensor);
        GreenFun(:,:,idx_sensor) = GreenFunGen(dGrid, dGrid, r_sensor,...
            2*Rd/Nd, 2*pi*fq/sound_speed);
    end
    GreenFun = permute(GreenFun, [3, 1, 2]);
    P_curr = reshape(P(:, fq_no), [num_sensor_points, 1, 1]);
    Phi = Phi - 1i*2*pi*fq * reshape(sum(...
        arrayfun(@(x,y) x.*y, GreenFun, P_curr), 1), [Nd, Nd]);
    if strcmp(PlotMode, 'on')
        imagesc(kgrid_recon.y_vec*1e3, kgrid_recon.x_vec*1e3, log(abs((Phi))));
		chndl = colorbar; title(['TBP, index of Fq = ', num2str(fq_no)]);
        colormap('hot'); ylabel('x-axis [mm]'); xlabel('y-axis [mm]');
		title(chndl, 'log'); set(gca,'fontsize', 16); axis square;
        drawnow
        F(fq_no) = getframe;
    end
end
if strcmp(PlotMode, 'on')
    movie2avi(F, 'TBP.avi', 'compression', 'None');
end
Phi = real(Phi')/num_sensor_points;
Phi = Phi/max(max(abs(Phi)));
if strcmp(GPUmode, 'on')
    Phi = gather(Phi);
end
end

function outMat = GreenFunGen(x_grid, y_grid, source, h, k)
% GreenFunGen
% Description: Calculate the free-space Green's function in 2D with
% smoothness on a singular point.
% Reference: JC Aguilar & Y Chen (Computers and Mathematics with
% Applications, 2014)
% Input:
%   x_grid, y_grid: grids on x-axis and y-axis
%   source: coordinate of source
%   h: radius of the covering range of a singular point
%   k: wave number
% Output:
%   outMat: the free-space Green's function
% Author: Hongxiang Lin
% Version: 1.0
% Date: 2014-12-09
[X,Y] = meshgrid(x_grid-source(1), y_grid-source(2));

c1 = -1.3105329259115095;
beta1const = -1i/4+1/2/pi*(log(h*k/2)-psi(1)+c1);

distMat = sqrt(X.^2+Y.^2).* (X.^2+Y.^2 > h^2)+ h*(X.^2+Y.^2 <= h^2);
outMat = 1i/4*((besselj(0, k*distMat)+1i*bessely(0, k*distMat)).* ...
    (X.^2+Y.^2 > h^2)+ beta1const*(X.^2+Y.^2 <= h^2));
end