%% Example of how to read the cell data outut from WavePart
%
% This is an example code to read the output parameters from the wavepart code.
% The various parameters for each partition and the bulk parameter are output
% as cells. For example the peak and mean wave frequencies for each partition
% as well as for the wahole spectrum are output in the f cell. The same applies
% for wave heights and wave directions. Below is a script that loads the data into 
% abd creates time-series. Also an example on how to plot teh frequencies as a 
% time-series is provided. The user shoudl adjust this script to fit her/his requirements
%
%
%% Load the datafile 
%
load ..\data\waveparamsoutput.mat  % replace with your datafile that stores the results of wavepartition
%  
%% Create time series from each data containing cell
%
for i = 1:length(f)     % Do loop for each time-step
    %
    % Mean (fm) and Peak (fp) frequency of each partition (Hz)
    %
    ff          = f{i};        % Extract cell content for this time-step i
    [m,n]       = size(ff);    % n-1 is the no of partitions
    fm(i,1:n-1) = ff(1,1:n-1); % take all fmean values for all partitions
    fp(i,1:n-1) = ff(2,1:n-1); % take all fpeak values for all partitions
    %
    % Mean (Dm) and Peak (Dp) direction (degs) and Directional spead (sigma) for each partition.
    %
    dd          = D{i};
    Dm(i,1:n-1) = dd(1,1:n-1);
    Dp(i,1:n-1) = dd(2,1:n-1);
    si(i,1:n-1) = dd(3,1:n-1);
    %
    % RMS (Hrms), significant (Hsig) wave height (m) and significant slope (Hsig/L)
    %
    hh            = H{i};
    Hrms(i,1:n-1) = hh(1,1:n-1);
    Hsig(i,1:n-1) = hh(2,1:n-1);
    psi(i,1:n-1)  = hh(3,1:n-1);
    %
    % Bulk parameters (from full (no-partitioned) specta)
    %
    fbm(i)   = ff(1,n); 
    fbp(i)   = ff(2,n); 
    Dbm(i)   = dd(1,n);
    Dbp(i)   = dd(2,n);
    sib(i)   = dd(3,n);
    Hbrms(i) = hh(1,n);
    Hbsig(i) = hh(2,n);
    psib(i)  = hh(3,n);
end
%
% Replace zero values (i.e., no partition identified) with NaN
fm(fm==0)=NaN;                
fp(fp==0)=NaN;
Dm(Dm==0)=NaN;
Dp(Dp==0)=NaN;
si(si==0)=NaN;
Hrms(Hrms==0)=NaN;
Hsig(Hsig==0)=NaN;
psi(psi==0)=NaN;

%% Example of plotting
%
figure
plot(fm(:,1:2)) % Plot values for partitions 1 (swell) and 2 (wind wave)
hold on
plot(fp(:,1:2))
legend('fm - #1','fm - #2','fp - #1', 'fp - #2')
title('mean and peak frequency for partitions 1 and 2')
ylabel('Frequency (Hz)')