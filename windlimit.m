%% windlimit.m
function [Fw,Fwpk] = windlimit(f,Sf,n)
% [Fw,Fwpk] = windlimit(f,Sf,n)
%
% Function that estimates the lower frequency limit of the wind wave spectra
% directly from the 2-D wave spectrum when no wind estimate is available. 
% It assumes that the wind spectra follows an f^-n (n=3 or 4, default = 4) 
% roll-off pattern. 
%
%% Inputs
%  f   - frequency array of spectrum (Hz)
%  Sf  - Spectral energy
%  n   - roll-off energy exponent (default n=4)
%
%% Output
%  Fw   - lowest frequency of wind-induced waves
%  Fwpk - frequency of wind wave peak energy 
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
%% Find Fw, lowest frequency of wind-induced waves
if nargin<3
    n=4;
end
% Searches for wind peaks within [fmin fmax] for more accurate estimates of
% cutoff, these do not have to be exact, but should be above the swell
% frequency ranges
f1 = 0.35;                      % Lower freq of initial wind range
f2 = 0.70;                      % Upper freq of initial wind range
df = f(2)-f(1);                 % Frequency resolution
S  = Sf/(sum(Sf)*df);           % Normalise the spectrum by total variance
iw = find(f>=f1 & f<=f2);       % Identify freq. range representing wind
yni = log10(S(iw).*f(iw).^n);   % estimate S(f)*f^n
ain = polyfit(f(iw),yni,1);     % fit a line on initial wind range
Sn  = polyval(ain,f);           % fitted line
S10 = log10(S.*f.^n);           % Log10 of fitted line
Dev = (S10(2:end)-Sn(2:end));   % Absolute value of deviation from S(f)*f^n
Varw= std(Dev(iw));             % std of deviation in the initial wind range
Limw= mean(Dev(iw));            % Mean value of deviation over the initial wind range
LIM = Limw-3*Varw;              % difference between mean and 3*std
X = cumsum(Dev-LIM);            % cumsum of deviation - LUM
if1 = find(f>f2,1);             % index of upper initial wind range, f2
[~,ik]=min(X(1:if1));           % look for wind minimum up to f2
if isempty(ik)
    ik=1;                       % if no min, return lowest frequency index
end
Fw  = f(ik);                    % Fw, lowest frequency of wind-induced waves
%
%% Find Fwpd frequency of wind wave peak energy - in frequency range (f>Fw)
i = f > Fw;
S1 = S.*i;
[~,i] = max(S1);
Fwpk = f(i);                    % frequency at peak wind wave energy
end