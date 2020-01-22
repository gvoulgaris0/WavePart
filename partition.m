%% partition.m
function [AA,Ef]=partition(freq,dir,E,wfc,fw,sw)
%%
% [AA,Ef]=partition(freq,dir,E,[wfc],[fw],[sw])
%
% Function to partition a directional wave spectrum in to different
% components. The partitions are given in energy descenting order.
% Partitions identified by the watershed function are modified following
% mostly the method described in:
%
% Hanson and Phillips (2001) Automated Analysis of Ocean Surface Directional
% Wave Spectra. Journal of Oceanic and Atmospheric Technology, 18, 278-293.
%
%
%% Inputs
%  freq  - Frequency array of directional spectrum (Hz)
%  dir   - Directional array of spectrum (degs)
%          NOTE: the same direction should not appear twice (i.e. no 0 & 360o
%          or +/-180o)
%  E     - Directional wave spectral density m2/Hz/deg
%  wfc   - (optional) (=1) only keep wind partitions within wind band
%  fw    - (optional) wind frequency lower limit (Hz), ex. 0.8*fpeak of wind
%  sw    - (optional) switch to plot (only if sw=999)
%
%% Outputs
%  AA    - Matrix indicating the partition each E(f,theta)value belongs to
%          =0 is the noise, =1 is the wind partition, >=2 are the swell partitions in
%          descending order of energy, an example with 2 swell partitions
%          Partition # - Partition type
%                    0 - Noise
%                    1 - Wind waves
%                    2 - First (more energetic) swell partition
%                    3 - Second (less energetic) swell partition 
%  Ef    - Smoothed energy matrix used for partition calculations,
%           
%% Uses
%  watreshed_ww3.m     - computes a matrix identifying the watershed
%                        regions of the input matrix (available as a mex
%                        file too for increased computational speed).
%  filterDirWavespec.m - smooths the measured spectrum (E)and creates Ef
%  peakspread.m        - peak spreading (df2) calculation as in Hanson and Phillips (2001)
%  valley_min.m        - lowest valley between partitions as in Hanson and Phillips (2001)
%  polarPcolor.m       - pcolor in polar coordinates developed by E. Cheynet (2019).
%                         (https://www.mathworks.com/matlabcentral/fileexchange/49040-pcolor-in-polar-coordinates)
%                         MATLAB Central File Exchange. Retrieved March 16, 2019.
%
%% Updates
% 
%  1/22/2020:  Added check for minimum energy required for at least one
%              partition to be generated. This avoids flat line spectra.
%
%  1/18/2020:  Added comment about data input should not contain the same
%              direction twice and other minor changes. Also lines
%              ee=zeros(N,N) and dd=zeros(N,N) were added.
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
if nargin<6 || isempty(sw)
    sw = 10;
end
if nargin<4 || isempty(wfc)
    sw  = 10;
    wfc = 1;
end
%% Set parameters for analysis
navg       = 3; % window size, ex. 3 uses averaging window of [3x3] twice, 5 uses [5x5] twice
% wind parameters
wind_width = 90; % width of wind band in degrees, default 120 degrees
windminf   = 0.12; % merge wind partitions above this frequency (Hz, default 0.15)
% swell parameters
swell_hs_lim = 0.2; % Hsig (m) minimum of swell partitions <---- most useful parameter to change
d_lim      = 90; % maximum number of degree separation for merging partitions
df         = freq(2)-freq(1);
dth        = dir(2)-dir(1);
minSqDist  = (6*df)^2;
% More swell parameters from Hanson and Phillips 2001 (See Table 1)
% to define the noise limit: ep > A/(fp^4+B)
%        N.Pac GofM GofM M.Pac S.Pac.
% A    =[60.00,2.00,2.00, 6.00,40.00]*1e-6;
% B    =[20.00,3.00,3.00,30.00,20.00]*1e-3;
% kappa=[ 0.40,0.40,0.40, 0.40, 0.50];
% z    =[ 0.65,0.65,0.65, 0.65, 0.75];
%
% very sensitive (up to 10 partitions)
% A = 2e-5;    %
% B = 0.05;    %
%
% less sensitive (2-4 partitions)
A = 4e-5;    % 12e-5;    %
B = 0.04;    %
kappa = 0.4; % Df2 < kappa*df2
% minimum between peaks EM > z*El where El is the smaller Ep of the two partitions
z = 0.4;     % 0.65; 
%
% optional input windminf
if nargin> 4
    windminf = fw; % input windminf (fw)
end
%% STEP 0: Check spectrum quality for partioning - 1/22/2020 Update
% A requirement is set that the spectrum has the minimum reuired energy for a single partition
%
minE_allowed = A/(max(freq).^4+B); % min energy for peak at highest frequency of data
if sum(E(:) - min(E(:))) < minE_allowed 
    error('spectrum has not enough energy for partitioning')
end

%% STEP 1: Filter measured spectra using double convolution
%
Ef = filterDirWavespec(E,2,navg);
%
%% STEP 2: Identify all partitions possible
%
% AA = watershed_ww3(Ef);        % WW3 watershed algorithm function (using mex is faster)
AA      = watershed_ww3_mex(Ef); % WW3 watershed algorithm (suggested mex)
[m,n]   = size(AA);
N       = double(max(max(AA)));   % number of partitions identified
% calculate partition parameters
fp = ones(N,1);
Dp = fp;
Ep = fp;
for i=1:N
    Mask1   = AA == i;
    [~,k]   = max(Ef(:).*Mask1(:));
    [I1,J1] = ind2sub(size(AA),k);
    fp(i)   = freq(I1);   % Peak frequency of partition i
    Dp(i)   = dir(J1);    % Peak direction of partition i
    Ep(i)   = Ef(I1,J1);  % Peak energy of partition i
end
%
%%   STEP 3: Define wind partition (wfc option) and merge all partitions within wind region
%
Ep(fp<windminf) = 0;    % ignore partitions with fp < the wind frequency mininimum (windminf)
if sum(Ep) < 0          % if there are no wind partitions throw error
    error('no wind partitions, try reducing input fw (ex. 0.1 Hz)')
end
[~,wN] = max(Ep);  % look for largest wind partition (fp > windminf)
fpw    = fp(wN);   % peak frequency of wind partition
dpw    = Dp(wN);   % peak direction of wind partition
% merge frequencies above this that are within cos(d-dtheta),
delta_th                = dir - dpw; % d-dtheta
delta_th(delta_th>180)  = delta_th(delta_th>180) - 360;
delta_th(delta_th<-180) = delta_th(delta_th<-180) + 360;
id                      = abs(delta_th) < wind_width; % locations within wind_width (default 90)
fc                      = windminf./cosd(delta_th*90/wind_width); 
fc(~id)                 = nan;
for i=1:N                            % merge partitions inside parabola
    if fp(i) == fpw
        continue
    else
        dj = find(dir == Dp(i),1);
        if fp(i) > fc(dj)            % if peak energy of partition is within wind parabola
            [parti]   = find(AA==i); % other partition number
            AA(parti) = wN;          % wind partition number
        end
    end
end
if wfc == 1                  % (wfc option) merge all points inside parabola and remove points outside parabola of the wind partition
    BB  = zeros(m,n);        % temp noise matrix
    for i = 1:length(dir)
        fcd = fc(i);
        if ~isnan(fcd)
            j = freq > fcd;  % keep points inside
            AA(j,i) = wN;
            k = freq < fcd;  % for the wind partition, remove pts outside
        else
            k = freq > -1;
        end
        BB(k,i) = 255;
    end
    [parti]   = find(AA==wN & BB == 255); % wind partition number noise values
    AA(parti) = 0; % set noise to 0
end
% Re-number wind partition to partition 1, all other partitions >2 are swell, noise = 0
BB    = AA;
ic = 1;
for i  = 1:max(max(AA))
    in = length(find(AA==i));
    if in~=0
        if i == wN % wind wind partition number =1
            [j]   = find(AA==i);
            BB(j) = 1;
        else
            ic    = ic+1;          % partition counting
            [j]   = find(AA==i);
            BB(j) = ic;
        end
    end
end
AA  = BB;
%
%%   STEP 4: Merge remaining swell partitions that are too close together (df2 < k*Df2 or minSqDist < (6*df)^2) and less than d_lim apart in direction
%
N  = max(max(AA));
fp = zeros(N,1);
Dp = fp;
df2= fp;
for i=2:N % ignore noise partition(=0) and wind partition(=1)
    Mask1   = AA == i;
    [~,k]   = max(Ef(:).*Mask1(:));
    [I1,J1] = ind2sub(size(AA),k);
    fp(i)   = freq(I1);   % Peak frequency of partition i
    Dp(i)   = dir(J1);    % Peak direction of partition i
    df2(i)  = peakspread(freq,dir,Ef,Mask1); % normalized spreading of partition i
end
% Find squared distances of partitions (SqDist)
fx      = fp.*cosd(Dp);
fy      = fp.*sind(Dp);
SqDist  = zeros(N,N);
Sqsprd  = zeros(N,N);
for i=2:N
    Df2         = (fx(i)-fx).^2+(fy(i)-fy).^2;
    SqDist(:,i) = Df2;
    Sqsprd(:,i) = ones(size(Df2))*df2(i);
end
% Merge partitions below a threshold (minSqDist)
% or if SqDist < kappa*Sqsprd   Hanson and Phillips 2001 Df2 < k*df2
for i=2:N
    L = find(SqDist(:,i) < minSqDist | ...
        SqDist(:,i) < kappa*Sqsprd(:,i));
    for d = 1:length(L)
        if L(d)~=i && L(d) > 1 % merge two partitions
            ddth = Dp(i) - Dp(L(d));
            ddth = abs(ddth);
            if ddth > 180
                ddth = ddth - 360;
            end
            if abs(ddth) < d_lim
                [j]    = find(AA==L(d));
                AA(j)  = i;
            end
        end
    end
end
% Re-number remaining swell partitions
BB = AA;
ic = 1;
for i  = 2:max(max(AA))
    in = length(find(AA==i));
    if in~=0
        ic    = ic+1;          % partition counting
        [j]   = find(AA==i);
        BB(j) = ic;
    end
end
AA = BB;
%
%% STEP 5: Select only partitions above noise level and partitions below 0.6 Hz (e <= A/(fp^4+B))
%
N = max(max(AA)); % number of partitions after merging
% find energy included in each partition.
P = zeros(N,1);
fp = P;
for i = 2:N                                % for each partition
    Mask1 = AA == i;
    P(i)   = sum(sum(Ef.*Mask1))*df*dth;   % Energy of partition class i
    [~,k]  = max(Ef(:).*Mask1(:));
    [I1,~] = ind2sub(size(AA),k);
    fp(i)  = freq(I1);                     % Peak frequency of partition i
end
jo    = 1:N;
mv    = A./(fp.^4+B);
jn    = find(P(:)>=mv(:) & fp(:) < 0.6 );  % Find partition classes containing > minEpart energy
op    = find(~ismember(jo,jn));            % Identify the partition index that passes the minVar criterion
for i  = 1:length(op)                      % Assign all low enery partitions (<minVar) to partition 0
    if op(i) > 1                           % ignore wind partition
        [j]    = find(AA==op(i));
        AA(j)= 0;                          % set to noise (=0)
    end
end
N   = length(jn);      % Number of remaining swell partitions
BB    = AA;
for i = 1:N            % Re-number remaining, energetic partitions starting at 2
    [j] = find(AA==jn(i));
    BB(j) = i + 1;
end
AA    = BB;
%
%%   STEP 6: Merge remaining swell partitions that do not have a valley between them and are separated
%    by less than d_lim (def 90) degrees and f_lim in Hz, and are next to
%    each other
%
N = max(max(AA));  % number of partitions after merging
fp = zeros(N,1);
Dp = fp;
Ep = fp;
for i=2:N          % ignoring noise partition(=0) and wind partition(=1)
    Mask1   = AA == i;
    [~,k]   = max(Ef(:).*Mask1(:));
    [I1,J1] = ind2sub(size(AA),k);
    fp(i)   = freq(I1);   % Peak frequency of partition i
    Dp(i)   = dir(J1);    % Peak direction of partition i
    Ep(i)   = Ef(I1,J1);  % Peak energy of partition i
end
vmin = zeros(N,N); % minimum energy on the line connecting two partitions peak energies
dd = zeros(N,N);   % added 1/20/2020
ee = zeros(N,N);   % added 1/20/2020
for i=2:N
    for ii = 2:N
        if i~=ii
            f1 = fp(i);
            f2 = fp(ii);
            d1 = Dp(i);
            d2 = Dp(ii);
            ee(ii,i) = min([Ep(i) Ep(ii)]);
            dd1 = abs(d2 - d1);
            if dd1 > 180
                dd1 = dd1 - 360;
            end
            dd(ii,i) = dd1;
            % minimum energy on the line connecting two partitions peak energies
            vmin(ii,i) = valley_min(Ef,freq,dir,d1,d2,f1,f2,AA);
        end
    end
end
% Merge partitions that do not have a valley between them (lower than z*ee,
% where ee is the lower peak energy of the two partitions, z is an
% empirical constant) and are separated by less than d_lim degrees
for i=2:N
    L = find((vmin(:,i) > z*ee(:,i)  & abs(dd(:,i)) < d_lim));
    for d = 1:length(L)
        if L(d)~=i && L(d) > 1 % merge different swell partitions
            [j]    = find(AA==L(d));
            AA(j)  = i;
        end
    end
end
% Re-number remaining partitions
N = max(max(AA));
BB    = AA;
ic = 1;
for i  = 2:N
    in = length(find(AA==i));
    if in~=0
        ic    = ic+1;                % partition counting
        [j]   = find(AA==i);
        BB(j) = ic;
    end
end
AA  = BB;
%
%% STEP 7 keep only swell partitions with significant waveheight above swell_hs_lim
N     = max(max(AA));                % number of partitions
sumSw   = zeros(N,1);
Hsig    = fp;
for i=2:N                            % only swell partitions
    Mask1   = AA == i;
    Sw      = E.*Mask1;
    sumSw(i)= sum(sum(Sw))*df*dth;   % Energy of partition class i
    Hsig(i)    = 4*sqrt(sumSw(i));   % Hsig wave height (m)
end
% keep only swell partitions with (default Hsig > 0.2 meters)
for i = 2:N
    if Hsig(i) < swell_hs_lim
        Mask1   = AA == i;
        AA(Mask1) = 0;               % set to noise
    end
end
% Re-number remaining partitions
N = max(max(AA));
BB    = AA;
ic = 1;
for i  = 2:N
    in = length(find(AA==i));
    if in~=0
        ic    = ic+1;                % partition counting
        [j]   = find(AA==i);
        BB(j) = ic;
    end
end
AA  = BB;
%
%% STEP 8: Re-order partitions in terms of energy level in descending order 
% Partition # ; Partition type
%    0        ;   Noise
%    1        ;   Wind
%    2        ;   largest (Hsig) swell partition
%    3        ;   2nd largest Swell partition
%
N = max(max(AA));
sumSw = zeros(N-1,1);
for i=2:N                              % for each partition
    Mask1 = AA == i;
    Sw      = E.*Mask1;
    sumSw(i-1)= sum(sum(Sw))*df*dth;   % Energy of partition class i
end
%
[~,Js] = sort(sumSw,'descend');        % Sort according to energy
BB         = AA;
for i=1:N-1
    [j]   = find(AA==Js(i)+1);
    BB(j) = i+1;
end
AA=BB;
%
%% Plot option
if sw == 999
    figure
    subplot(121)
    surf(freq',dir',Ef',AA')
    title([num2str(Nw2) ' partitions'])
    subplot(122)
    [~,c] = polarPcolor(freq',[dir ; dir(1)]',[AA  AA(:,1)]');
    c.Ticks = 1:Nw2;                  % np ticks
    cm = colormap;                    % 64 default colors
    cm = cm(1:64/Nw2:64,:);           % reduce to np colors
    colormap(cm)
end
%
end % end function