%% waveparamspart.m
function [f,D,Ee,H] = waveparamspart(E,freq,dir,AA,h)
%%  
%  [f,D,Ee,H] = waveparamspart(E,freq,dir,AA,[h])
%
%  Function to calculate the wave parameters for each partition of the
%  spectrum. The peak values for each partition are estimated as well as
%  the mean and other parameters as described in Hanson and Phillips (2001)
%  (see Appendix therein).
%
%% Inputs
%  E    = Directional Wave energy density (m2/Hz/deg)
%  freq = frequency array of spectral estimates (Hz)
%  dir  = direction array of spectral estimates (degs)
%  AA   = Matrix same dimensions as E indentifying partition element by a number.
%  h    = water depth [Optional]
%
%% Outputs
%  f(2,n) = [fm fp] Mean and Peak frequency of each partition n (Hz)
%  D(3,n) = [Dm Dp sigma] Mean and Peak direction (degs) and Directional spead of each partition n
%  Ee(2,n)= [Et Ep] Total and Peak energy of each partition n (m2 and m2/Hz/degs)
%  H(3,n) = [Hrms Hsig psi] rms and significant wave height (m) and significant slope
%
%% Uses 
%  dispersion.m  - function to solve the dispersion equation for shallow waters
%
%% Authors
%  Douglas Cahl and George Voulgaris
%  School of the Earth, Ocean and Environment
%  University of South Carolina, Columbia, SC, USA
%
%% Updates
%  01/21/2020 - Corrections in the formulas estimating Hrms and Dm
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
if nargin<5 || isempty(h)
    h=1000;   % Deep water conditions
end
g   = 9.81;
Nw1 = max(max(AA)); % total number of partitions
f   = ones(2,Nw1);
D   = ones(3,Nw1);
Ee  = f;
H   = D;
df  = freq(2)-freq(1);
dth = dir(2)-dir(1);
[nf,~] = size(E);   % size of the 2-d matrix directional spectra array  
                    % nf = no of freq. bins, nd (~)= no of direction bins
the  = repmat(dir(:)',nf,1);
%
for i = 1:Nw1                                 % for each partition i
    Mask1   = AA == i;                        
    Epart   = E.*Mask1;                       % Si(f,theta)
    [~,k]   = max(Epart(:));
    [I1,J1] = ind2sub(size(AA),k);
    fp      = freq(I1);                       % Peak frequency for Si
    Dp      = dir(J1);                        % Peak direction for Si
    Ep      = E(I1,J1);                       % Peak energy level for Si
    Et      = sum(sum(Epart))*df*dth;         % \int Si(f,theta)
    Hrms    = 2*sqrt(2*Et);                   % Hrms wave height (m) - 1/21/2020
    Hsig    = 4*sqrt(Et);                     % Hsig wave height (m)
    fm      = sum(sum(Epart,2)./freq)*df*dth; % Mean freq. of Si
    f(1:2,i)= [fm fp];
    sino    = sum(sum(Epart.*sind(the)))/Et;  % normalised
    coso    = sum(sum(Epart.*cosd(the)))/Et;  %
    Dm      = atan2d(sino,coso);              % Mean Direction of Si - 1/21/2020
    epsilon = sqrt(1-(sino^2+coso^2));
    sigma   = (1+0.1547*epsilon^3)/sin(epsilon);% Directional spread of Si
    D(1:3,i)= [Dm,Dp,sigma];
    Ee(1:2,i)= [Et,Ep];
    kh      = dispersion((2*pi*fp).^2*h/g);
    Lp      = 2*pi/(kh/h);
    psi     = Hsig/Lp;                        % Significant slope (Huang 1986)
    H(1:3,i)= [Hrms,Hsig,psi];
end
end