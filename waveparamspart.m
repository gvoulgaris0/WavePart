%% waveparamspart.m
function [f,D,Ee,H,Er,Np] = waveparamspart(E,freq,dir,AA,h)
%%  
%  [f,D,Ee,H,Er,Np] = waveparamspart(E,freq,dir,AA,[h])
%
%  Function to calculate the wave parameters for each partition of the
%  spectrum plus the total bulk parameters using the full, unpartitioned spectrum. 
%  The peak and mean values for each partition are estimated as well as
%  other parameters as described in Hanson and Phillips (2001),(see Appendix therein).
%
%% Inputs
%  E    = Directional Wave energy density (m2/Hz/deg)
%  freq = frequency array of spectral estimates (Hz)
%  dir  = direction array of spectral estimates (degs)
%  AA   = Matrix same dimensions as E indentifying partition element by a number.
%  h    = water depth [Optional]
%
%% Outputs
%  f (1:2,1:Np+1) = [fm fp] Mean and Peak frequency of each partition (Hz)
%  D (1:3,1:Np+1) = [Dm Dp sigma] Mean and Peak direction (degs) and Directional spead of each partition n
%  Ee(1:2,1:Np+1) = [Et Ep] Total and Peak energy of each partition n (m2 and m2/Hz/degs)
%  H (1:3,1:Np+1) = [Hrms Hsig psi] rms and significant wave height (m) and significant slope
%  Er             = Diagnostic error (in %) showing the energy not represented by
%                   the partitions. Er= 100*( Etotal-sum(Epartitions) ) / Etotal.
%  Np             = Number of partitions
% 
%  note:  the bulk parameter values are reported at location Np+1, after
%  the partitioned parameters.
%
%% Uses 
%  dispersion.m  - function to solve the dispersion equation for shallow waters
%
%% Updates
%  10/13/2019    -  The diagnostic Er and bulk parameter estimations were
%                   added. Also the mean period estimation was wrong and it
%                   has been corrected. 
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
if nargin<5 || isempty(h)
    h=1000;   % Deep water conditions
end
g   = 9.81;
Np  = max(max(AA)); % total number of partitions
f   = ones(2,Np+1);
D   = ones(3,Np+1);
Ee  = f;
H   = D;
df  = freq(2)-freq(1);
dth = dir(2)-dir(1);
[nf,~] = size(E);   % size of the 2-d matrix directional spectra array  
                    % nf = no of freq. bins, nd (~)= no of direction bins
the  = repmat(dir(:)',nf,1);
%
for i = 1:Np+1                               % for each partition i, plus total
    if i<=Np
        Mask1 = AA == i;
        Epart = E.*Mask1;                     % Partition i (Si(f,theta))
    else
        Epart = E;                            % Total (S(f,theta))
    end
    [~,k]   = max(Epart(:));
    [I1,J1] = ind2sub(size(Epart),k);
    fp      = freq(I1);                       % Peak frequency for Si
    Dp      = dir(J1);                        % Peak direction for Si
    Ep      = E(I1,J1);                       % Peak energy level for Si
    Et      = sum(sum(Epart))*df*dth;         % \int Si(f,theta)
    Hrms    = sqrt(2*Et);                     % Hrms wave height (m)
    Hsig    = 4*sqrt(Et);                     % Hsig wave height (m)
    fm      = sum(sum(Epart,2).*freq)/...     % Mean freq. of Si (m1/m0)
              sum(sum(Epart));
    f(1:2,i)= [fm fp];
    sino    = sum(sum(Epart.*sind(the)))/Et;  % normalised
    coso    = sum(sum(Epart.*cosd(the)))/Et;  %
    Dm      = atan(sino/coso);                % Mean Direction of Si
    epsilon = sqrt(1-(sino^2+coso^2));
    sigma   = (1+0.1547*epsilon^3)/sin(epsilon);% Directional spread of Si
    D(1:3,i)= [Dm,Dp,sigma];
    Ee(1:2,i)= [Et,Ep];
    kh      = dispersion((2*pi*fp).^2*h/g);
    Lp      = 2*pi/(kh/h);
    psi     = Hsig/Lp;                        % Significant slope (Huang 1986)
    H(1:3,i)= [Hrms,Hsig,psi];
end
Er=100*(Ee(1,Np+1)-sum(Ee(1,1:Np),2))/Ee(1,Np+1); % Error in total energy partitioned (in %)
end