%% valley_min.m
function vmin = valley_min(E,freq,dir,d1,d2,f1,f2,AA)
%%
%  Function used to find the lowest point (valley) between two local peaks
%  in the frequency directional spectrum.
%
%  Inputs
%    E    = 2-D matrix (dir - freq.) to be searched
%    freq = frequency array for E 
%    dir  = direction array for E
%    d1   = direction of partition 1
%    d2   = direction of partition 2
%    f1   = frequency of partition 1
%    f2   = frequency of partition 2
%    [AA] = partitions, noise = 0, optional
%
%  Output
%    vmin = the minimum value of the line connecting (d1,f1) to
%    (d2,f2), returns Nan if passes through noise (AA=0) (optional input)
%    or partitions different than partition 1 or 2
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

% indexes of partition peaks
i1 = find(freq==f1);
j1 = find(dir==d1);
i2 = find(freq==f2);
j2 = find(dir==d2);
    
% dir runs from -180 to 180
if abs(d1 - d2) > 180 % we need to wrap around instead
    % reference from point d2
    j2 = find(dir==d2);
    if d2 > d1 % wrap points from right to left
        % [-180 ... d1 ... d2 ... 180] -> [d2-360 ... -180 ... d1 ... ]
        if j2 > 1
            d2  = d2 - 360;
            E   = [E(:,j2:end) E(:,1:j2-1)];
            dir = [dir(j2:end)-360; dir(1:j2-1)];
        else
            E  = E(:,j2:end);
            dir = dir - 360;
        end
    else % wrap points from left to right
        % [-180 ... d2 ... d1 ... 180] -> [... d1 ... 180 ... d2 + 360]
        if j2 > 1
            d2  = d2 + 360;
            E   = [E(:,j2:end) E(:,1:j2-1)];
            dir = [dir(j2:end); dir(1:j2-1)+360];
        else
            E  = E(:,j2:end);
            dir = dir + 360;
        end
    end
end
ln_dir = linspace(d1,d2);               % direction pts
m_val  = (f2-f1)/(d2-d1);               % slope
ln_f   = f1 + m_val*(ln_dir-d1);        % line between (d1,f1) and (d2,f2)
[f,th] = meshgrid(freq(:),dir(:));      % meshgrid 
ln_val = interp2(f,th,E',ln_f,ln_dir);  % values on line
vmin = min(ln_val);                     % minimum value on line
if nargin > 7 && ~isnan(m_val) % check if crosses noise partition (AA=0)
    P1 = AA(i1,j1);
    P2 = AA(i2,j2);
    for i = 1:length(ln_dir)
        d = ln_dir(i);
        [dcheck,di] = min(abs(dir - d));
        if dcheck > 20
            error('check valley_min.m')
        end
        f = ln_f(i);
        [~,fi] = min(abs(freq - f));
        atest = AA(fi,di);
        if atest == 0 || (atest ~= P1 && atest ~= P2)
            vmin = nan;
            break
        end
    end  
end
end

