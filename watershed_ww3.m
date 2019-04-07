%% watershed_ww3.m
function AA = watershed_ww3(E)
%%
% AA = watershed_ww3(E)
%
% Function that computes a matrix identifying the watershed regions of the
% input matrix E. It is similar to the Matlab(r) buil in function watershed.m
% The elements of AA are positive integers >=1 and each number labels a wateshed.
% The function uses 8-connected neighborhood points to identify the
% watersheds and allows for continuity between the boundaries.
%
% Input
%   E - 2D energy spectrum E(freq,Dir)
%
% Output
%   AA(N-freq,M-dir) - watershed partition number
%   Np               - Total number of partitions
%
%% Uses
%   nextmax.m - internal function
%
%% Authors
%  Douglas Cahl
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
AA     = zeros(size(E));    % peak matrix
[N,M]  = size(E);
pi     = nan(length(E(:)),1);
pj     = nan(length(E(:)),1);
pE     = nan(length(E(:)),1);
% ----------- watershed, doesn't deal with flat surfaces
Np = 0;
for i = 1:N     % for each frequency
    for j = 1:M % for each direction
        i0 = i;
        j0 = j;
        while true
            % find surrounding highest point
            [i1,j1,flat] = nextmax(E,i0,j0);
            if i0 == i1 && j0 == j1 && flat == 0  % if this is the peak
                if Np == 0                        % if first peak
                    Np = Np + 1;                  % part = 1
                    AA(i,j) = Np;                 % part matrix
                    pi(Np) = i0;                  % i value of partition peak
                    pj(Np) = j0;                  % j value of partition peak
                    pE(Np) = E(i0,j0);            % value of partition peak
                else
                    partnew = 1;                  % Check if this is a new peak
                    for k = 1:Np                  % if it belongs to peak already found
                        if pi(k) == i0 && pj(k) == j0
                            AA(i,j) = k;
                            partnew = 0;          % this is not a new peak
                            break;
                        end
                    end
                    if partnew == 1               % only runs this if loops above do not find a peak
                        Np = Np + 1;              % make new parition
                        AA(i,j) = Np;
                        pi(Np) = i0;
                        pj(Np) = j0;
                        pE(Np) = E(i0,j0);
                    end
                end
                break
            else                                  % run again to highest surrounding point in pt0
                i0 = i1;
                j0 = j1;
            end
        end
    end
end
end
% subroutine
function [i1,j1,flat] = nextmax(E,i,j)
N       = size(E,1);
M       = size(E,2);
% find highest point surrounding 22 in the matrix below
% incorporates wrap around conditions at the boundaries
% 31 32 33
% 21 22 23
% 11 12 13
pt      = nan(3,3);
pt(2,2) = E(i,j);

% same freq (row 2)
if j == M
    pt(2,3) = E(i,1);
else
    pt(2,3) = E(i,j+1);
end
%
if j == 1
    pt(2,1) = E(i,end);
else
    pt(2,1) = E(i,j-1);
end
% upper freq (row 1)
if i < N
    pt(3,2) = E(i+1,j);
    if j == M
        pt(3,3) = E(i+1,1);
    else
        pt(3,3) = E(i+1,j+1);
    end
    if j == 1
        pt(3,1) = E(i+1,end);
    else
        pt(3,1) = E(i+1,j-1);
    end
end
% lower freq (row 3)
if i > 1
    pt(1,2) = E(i-1,j);
    if j == M
        pt(1,3) = E(i-1,1);
    else
        pt(1,3) = E(i-1,j+1);
    end
    if j == 1
        pt(1,1) = E(i-1,end);
    else
        pt(1,1) = E(i-1,j-1);
    end
end
%
ptflat = pt(~isnan(pt));
ptflat = unique(ptflat);
if length(ptflat) == 1 % flat surface all points are the same
    flat = 1;
    % randomly push point somewhere
    if i == N
        i1 = randi(2);
    elseif i == 1
        i1 = randi(2) + 1;
    else
        i1 = randi(3);
    end
    i1 = i + i1 - 2;
    j1 = randi(3);
    if j == M && j1 == 3
        j1 = 1;
    elseif j == 1 && j1 == 1
        j1 = M;
    else
        j1 = j + j1 - 2;
    end
else
    flat = 0;
    % next highest point
    [~, index] = max(pt(:));
    [i1, j1] = ind2sub(size(pt), index);
    i1 = i + i1 - 2;
    if j == M && j1 == 3
        j1 = 1;
    elseif j == 1 && j1 == 1
        j1 = M;
    else
        j1 = j + j1 - 2;
    end
end
end