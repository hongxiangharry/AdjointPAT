function [Phi] = BackProjection(P, t_array, Nd, Rd, sensor, sound_speed, GPUmode)
% BackProjection
% Description: calculate the initial pressure distribution for PAT by using
% the Back-Projection method. ROI is a square centered at the origin.
% Input:
%   P: received waveform data
%   t_array: time sequence for the received waveform data
%   Nd: grid points per side of ROI
%   Rd: half length of a side of ROI
%   sensor: sensor information defined by the k-Wave toolbox
%   sound_speed: sound speed of the background
%   GPUmode: on/off of GPU mode
% Output:
%   Phi: the BP imaging function
% 
% Author: Hongxiang Lin
% Affiliation: Dept. Mechanical Engineering, the University of Tokyo
% 
% Version: 1.0
% Date: 2018-06-20
% 
if nargin < 7
    GPUmode = 'off';
end
dGrid = linspace(-Rd, Rd, Nd);

if strcmp(GPUmode, 'off')
    Phi = zeros(Nd, Nd);
    for idx_sensor = 1:size(sensor.mask,2)
        r_sensor = sensor.mask(:, idx_sensor);
        GreenFun = WaveEqGreenFunGen(t_array, r_sensor, dGrid, sound_speed);
        GreenFun = GreenFun(2:length(t_array),:,:)-GreenFun(1:length(t_array)-1,:,:);
        % debug
        for idx_t = 1:length(t_array)-1
            temp = reshape(GreenFun(idx_t, :, :), [Nd, Nd]);
            Phi = Phi + temp .* P(idx_sensor, idx_t);
        end
    end
elseif strcmp(GPUmode, 'on')
    Phi = zeros(Nd, Nd, 'single', 'gpuArray');
    for idx_sensor = 1:size(sensor.mask, 2)
        r_sensor = sensor.mask(:, idx_sensor);
        GreenFun = WaveEqGreenFunGen(gpuArray(single(t_array)), ...
                                      gpuArray(single(r_sensor)),...
                                      gpuArray(single(dGrid)),...
                                      gpuArray(single(sound_speed)));
        GreenFun = GreenFun(2:length(t_array),:,:)-GreenFun(1:length(t_array)-1,:,:);                                
        P_curr = reshape(P(idx_sensor, 1:length(t_array)-1), [length(t_array)-1, 1, 1]);
        Phi = Phi + reshape(sum(arrayfun(@(x,y) x.*y, GreenFun, P_curr),...
        1), [Nd, Nd]);
    end
    Phi = gather(Phi);
end
Phi(isnan(Phi)) = 0; % remove the exceptional nan-value.
Phi = real(Phi/max(max(abs(Phi))));
end

function GreenFun = WaveEqGreenFunGen(t_array, r_sensor, dGrid, sound_speed)
% WaveEqGreenFunGen
% Description: a time-domain Green's function generator for wave equation
% in homogeneous medium.
% Input:
%   t_array: time sequence for the received waveform data
%   r_sensor: coordinate of the sensor
%   dGrid: grid of a side of ROI
%   sound_speed: sound speed of the background
% Output:
%   GreenFun: the output Green's function
%
% Author: Hongxiang Lin
% Affiliation: Dept. Mechanical Engineering, the University of Tokyo
% 
% Version: 1.1
% Date: 2018-06-19
% History: 
% 2018-06-19: distMat is a 3d distance matrix. In GPU computing, sqrt can
% only address the positive real number. So abs is required for distMat. 
%
[T, X, Y] = ndgrid(t_array, dGrid, dGrid);
distMat = T.^2-((X-r_sensor(1)).^2+(Y-r_sensor(2)).^2)/sound_speed^2;
GreenFun = (distMat > 0).*(1./(2*pi)./sqrt(abs(distMat)));
end