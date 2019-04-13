%% wildpoint.m
function y = wildpoint(x,nn)
% y = wildpoint(x,nn)
%
% Wildpoint removal function. 
% If the value of point in the data series "x" differs by a value greater than
% nn x std from the median of the surrounding 4 points (2 on either side of it)
% it is replaced by a NaN
%
%% Inputs:   
%  x  = data series
%  nn = defines threshold (how many standard deviations)
%
%% Output:
%  y  = wildpoint cleaned data series
%
%% Authors
%  George Voulgaris
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
n = length(x);                   % length of data series
y = x;                           % preallocate y
for i = 1:n
    i1 = max(1,i-2);             % 2 points before
    i2 = min(n,i+2);             % 2 points after
    y(i) = median(x(i1:i2));     % 5 pt window median
end
xdiff      = abs(x - y);         % difference to median
xstd       = std(x);             % std of window
i_wild     = xdiff > xstd*nn;    % difference exceeeds threshold (nn*std)
y(i_wild)  = nan;                % replace wildpoints with NaNs
end