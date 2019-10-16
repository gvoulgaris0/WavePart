%% readspectra.m
function [fw,dw] = readspectra(t,freq,dir,S,swplot)
%
% Function that estimates the lower frequency limit of the wind wave
% spectra directly from a timeseries of 2-D wave spectra. For each wave
% spectrum, windlimit.m calculates the cutoff frequency and direction.
% After this, wildpoint removal and smoothing is applied to the frequency
% and direction calculations. See windlimit.m for more information.
%
%% Inputs
%  t      - array of time values
%  freq   - frequency array of spectrum (Hz)
%  dir    - direction array of spectrum (degrees)
%  S      - timeseries of 2-D wave spectral energies S(freq,dir,time)
%  swplot - >0 show plot of analysis
%
%% Output
%  fw   - lowest frequency of wind-induced waves (Hz, filtered+smoothed)
%  dw   - direction at peak wind wave energy (degrees, math orientation)  
%
%% Uses
%  windlmit.m, wildpoint.m
%
%% Authors
%  Douglas Cahl and George Voulgaris
%  School of the Earth, Ocean and Environment
%  University of South Carolina, Columbia, SC, USA
%
%% Copyright 2019 Douglas Cahl, George Voulgaris
%
% This file is part of WavePart.
%
% WavePart is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%% Main Function
%% Run the analysis on each wave spectrum
N   = length(t);  % length of timeseries
Fw  = zeros(1,N); % wind wave frequency lower limit (Hz)
Dpk = zeros(1,N); % direction (degrees) at peak wind wave energy (f>Fw)
Sf  = zeros(length(freq),N); % frequency spectrum
Sth = zeros(length(dir),N);  % direction specturm
for i = 1:length(t)
        E    = S(:,:,i);
        Sf(:,i)=conv(sum(E,2),[1 1 1 1 1]/5,'same');  % 5pt smoothed frequency spectrum
        Sth(:,i)=conv(sum(E,1),[1 1 1 1 1]/5,'same'); % 5pt smoothed direction spectrum
        [Fw(i),Fwpk] = windlimit(freq,Sf(:,i),4);
        j = freq==Fwpk; 
        Sdir = E(j,:);      % direction spectrum at peak wind frequency
        [~,k] = max(Sdir); 
        Dpk(i) = dir(k);    % direction at peak frequency
end
%
%% Wildpoint removal and smoothing
Const = 1; % how many stds in the 5pt median sort for wildpints remove
mvpt  = 15; % moving mean of window this wide for final result
% frequency (fw)
fw=wildpoint(Fw,Const);
fw=wildpoint(fw,Const);
j = ~isnan(fw);
fw = interp1(t(j),fw(j),t); % interp over nans
fw = movmean(fw,mvpt); % moving average
fw = movmean(fw,mvpt); % moving average
% direction (dw)
D = Dpk;
D(~j) = nan;
Dx = cosd(D);
Dy = sind(D);
Dx = interp1(t(j),Dx(j),t); % interp over nans
Dy = interp1(t(j),Dy(j),t); % interp over nans
Dx = movmean(Dx,mvpt);
Dx = movmean(Dx,mvpt); % moving average
Dy = movmean(Dy,mvpt);
Dy = movmean(Dy,mvpt); % moving average
dw = atan2d(Dy,Dx);
%
%% Plot (if swplot>0)
if swplot > 0
    figure
    % fw
    subplot(221) 
    plot(t(j),Fw(j),'.')
    hold on
    plot(t(~j),Fw(~j),'x')
    plot(t,fw,'-k','linewidth',2)
    legend(['good ' num2str(sum(j)) ],...
        ['removed ' num2str(sum(~j))],...
        [num2str(mvpt) ' pt avg'])
    ylabel('freq (Hz')
    datetick('x','keeplimits')
    title('wind wave low freq (fw)')
    % dw
    subplot(222) 
    plot(t(j),Dpk(j),'.')
    hold on
    plot(t(~j),Dpk(~j),'x')
    plot(t,dw,'-k','linewidth',2)
    legend('good','removed',[num2str(mvpt) ' pt avg'])
    ylabel('deg (math)')
    datetick('x','keeplimits')
    title('wind dir at peak wind energy (dw)')
    % Sf
    subplot(223) 
    pcolor(t,freq,log10(Sf))
    shading interp
    hold on
    plot(t,fw,'w-')
    ylabel('wave freq (Hz)')
    datetick('x','keeplimits')
    lg = legend('Sf','fw');
    lg.Color     = [0.1 0.1 0.1]; % white text on legend
    lg.TextColor = [1.0 1.0 1.0]; % black legend background
    % Sth
    subplot(224)
    pcolor(t,dir,log10(Sth))
    shading interp
    hold on
    plot(t,dw,'w-')
    ylabel('wave dir deg (math)')
    datetick('x','keeplimits')
    lg = legend('Sth','dw');
    lg.Color     = [0.1 0.1 0.1]; % white text on legend
    lg.TextColor = [1.0 1.0 1.0]; % black legend background
end
end