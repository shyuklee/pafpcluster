load('sampledata.mat');
% yx: (npts, 2) coordinate of localization [pixel]
%     (*, 1): y(column) coordiate
%     (*, 2): x(row) coordiate
% frFirstON: (npts, 1) first fluorescence ON frame
% frLastON : (npts, 1) last fluorescence ON frame

rMax = 1; % maximum nearest neighbor distance to be clustered [pixel] 
tF = 3.2; % [minute] time within which all Dendra2 get photoactibvated with uniform photoactivation rate
         % Fermi time when Fermi photoactivation scheme is employed.
secPerFrame = 0.0535; % second per frame
tcInit = 1; % initial tc valule [frame]

hFun = @(nMol) (optmtaucsolve2_dendra(nMol, tF)/secPerFrame); % optimal tc function handle [frame] 
[nCluster, idCluster] = optimalsptcluster(yx, frFirstON, frLastON, rMax, tcInit, hFun);
