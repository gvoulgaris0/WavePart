%% peakspread.m
function df2 = peakspread(freq,theta,E,mask)
%%
% df2 = peakspread(freq,theta,E,mask)
%
% Function to find the normalized spreading of each partition 
% as in eqn (7) in Hanson and Phillips 2001
% 
%% Inputs:
%  freq         = frequency array (1:nf)
%  theta        = direction array (1:nd)
%  E(nf,nd)     = the 2-D directional spectrum of the partition (i.e, Ep=E*Mask).
%  mask(nf,nd)  = array of 0s and 1s where the partition =1
%
%% Output:
%  df2          = squared spreading parameter \delta f^2
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
Ep      = E.*mask;
[nf,nd] = size(E);  % size of the 2-d matrix directional spectra array  
                    % nf = no of freq. bins, nd = no of direction bins
f   = repmat(freq(:),1,nd);
th  = repmat(theta(:)',nf,1);
e   = sum(sum(Ep)); 
Fx1 = (1/e)*sum(sum(Ep.*f.*cosd(th)));
Fy1 = (1/e)*sum(sum(Ep.*f.*sind(th)));
Fx2 = (1/e)*sum(sum(Ep.*f.^2.*cosd(th).^2));
Fy2 = (1/e)*sum(sum(Ep.*f.^2.*sind(th).^2));
df2 = Fx2 - Fx1^2 + Fy2 - Fy1^2;
end