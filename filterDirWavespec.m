%%  filterDirWavespec.m
function Ew = filterDirWavespec(E,nc,n)
%% 
%  Function used to apply a 2-D moving averaging filter on a matrix
%  representing directional spectrum. The objective is to allow that the
%  two direction ends converge from 0 to 360 degs and vice versa similarly
%  to repeat in boundary conditions in numerical modeling.
%
%  Inputs
%    E  = matrix to be convoluted
%    nc = number of convolutions (default = 1)
%    n  = size of convolution window (default = 3)
%
%  Output
%    Ew = filtered matrix E of the same dimensions
%
%  E =
%  |------------------------------|
%  |  |-----------------------|
%  |  |  -3  -2  1  1  1  1 -20 -30
%  20 30  3   2  1  1  1  1  20  30
%         3   2  1  1  1  1  20  30
%         3   2  1  1  1  1  20  30
%         3   2  1  1  1  1  20  30 3 2
%         4   3  1  1  1  1  25  35 | |
%         |-------------------------| |
%             |-----------------------|
%
%   with n=5, becomes
%
%  E2 =
%  -20 -30 -3 -2  1  1  1  1 -20 -30 -3 -2
%  -20 -30 -3 -2  1  1  1  1 -20 -30 -3 -2
%  -20 -30 -3 -2  1  1  1  1 -20 -30 -3 -2
%   20  30  3  2  1  1  1  1  20  30  3  2
%   20  30  3  2  1  1  1  1  20  30  3  2
%   20  30  3  2  1  1  1  1  20  30  3  2
%   20  30  3  2  1  1  1  1  20  30  3  2
%   25  35  4  3  1  1  1  1  25  35  4  3
%   25  35  4  3  1  1  1  1  25  35  4  3
%   25  35  4  3  1  1  1  1  25  35  4  3
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
if nargin<2
    n=3; nc=1;
end
if nargin<3
    nc=1;
end
if mod(n,2)== 0   % Even number
    n=n+1;        % Make it odd
end
W=ones(n,n)/(n^2); % n x n normalized window to be used
%[nf,nd]=size(E);  % size of the 2-d matrix array representing directional spectra
% nf = no of freq. bins, nd = no of direction bins
mn = floor(n/2);   % size of extra columns / rows the window will create
% repeat the 1st and last frequency columns at the start and end
%
E1 = [repmat(E(1,:),mn,1); E; repmat(E(end,:),mn,1) ];
% place the last mn directions to the top (in reversed order) and the first to the end in reversed order
E2 = [E1(:,end-mn+1:end), E1, E1(:,1:mn)];
for i=1:nc
    E3 = conv2(E2,W,'same');
    E2 = E3;
end
% remove the extra columns / rows that had been placed so Ew has the same size as E
Ew = E3(mn+1:end-mn, mn+1:end-mn);
end