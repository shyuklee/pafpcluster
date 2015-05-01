%+
% NAME: 
% 	optimalsptcluster
%
%
% PURPOSE:
% 	Clusterize spatio-temporal localization data into disinct molecules using 
% 	iterative optimal tc algorithm 
%
%
% CALLING SEQUENCE:
% 	[nCluster, idCluster] = optimalsptcluster(posArr, t1Arr, t2Arr, rMax, tc, hFun); 
%
%
% INPUTS:
% 	posArr: M x 1 (1D), M x 2 (2D),  M x 3 (3D),... M x N (ND) 
% 		position coordinate array [pixel]
% 	t1Arr: M x 1 first time coordinate array [frame]
% 	t2Arr: M x 1 last time coordinate array [frame]
% 	rMax: Maximum allowed spacial separation [pixel]
% 	tc: Initial tc value [frame]
% 	hFun: Function handle of computing optimal tc for a given number of molecules
% 	
%
% OUTPUTS:
% 	nCluster : Number of cluster
% 	idCluster : Cluster indices 
%
%
% PROCEDURE:
% 	Number of molecules is determined by temporal clustering ('fastsptcluster') 
% 	for a given tc value. An optimal tc value is determined by a user provided function
% 	'hFun' for a given number of molecules. These two main steps are repeated until 
% 	reaching convergence. 
% 	
%
% NOTES:
% 	Iterative optimal-tc spatio-temporal clustering 
%
%
% REFERENCE:
% 1. S. Lee, J. Y. Shin, A. Lee and C. Bustamante, 
%    "Counting single photoactivatable fluorescent molecules by 
%     photoactivated localization microscopy (PALM),"
%     Proc. Natl. Acad. Sci. 109, 17436-17441 (2012).
% 
%
% MODIFICATION HISTORY:
% 01/10/2012: Created by S.H.Lee, UC Berkeley 
%
%  Copyright (c) 2012 Sang-Hyuk Lee 
%
%
% LICENSE:
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%-

function [nCluster, idCluster] = optimalsptcluster(posArr, t1Arr, t2Arr, rMax, tc, hFun);

TOL = 0.01; 	% Convergence criteria 
iterMax = 50;	% Maximum number of iteration

% Initial clustering
[nCluster, idCluster] = fastsptcluster(posArr, t1Arr, t2Arr, rMax, tc); 
tOld = hFun(nCluster);

% Iterate
for i=1:iterMax
   [nCluster, idCluster] = fastsptcluster(posArr, t1Arr, t2Arr, rMax, tOld); % Clustering 
   tNew = hFun(nCluster); % Optimal tc
   if abs(tNew-tOld)/tOld < TOL
      [nCluster, idCluster] = fastsptcluster(posArr, t1Arr, t2Arr, rMax, ceil(tOld));
      break;
   else
      tOld = tNew;
   end
end % end for

% Display the result
disp(sprintf('Iteration: %d,   Tau_c: %.1f,   nBurst: %d,   nCluster: %d \n', i, tNew, numel(t1Arr), nCluster));
