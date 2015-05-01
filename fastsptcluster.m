%+
% NAME: 
% 	fastsptcluster
%
% PURPOSE:
% 	Fast version of sptcluster. Clustering order may be different from sptcluster.
%
% CATEGORY:
%
% CALLING SEQUENCE:
% 	[nCluster, idCluster] = fastsptcluster(posArr, t1Arr, t2Arr, rMax, tMax); 
%
% INPUTS:
% 	posArr : M x 1 (1D), M x 2 (2D),  M x 3 (3D),... M x N (ND) 
% 		position coordinate array [pixel] 
% 	t1Arr : M x 1 first time coordinate array [frame]
% 	t2Arr : M x 1 last time coordinate array [frame] 
% 	rMax : Maximum allowed spacial separation [pixel]
% 	tMax : Maximum allowed temporal separation [frame]
% 	
%
% OUTPUTS:
% 	nCluster : Number of cluster
% 	idCluster : Cluster indices 
%
%
% PROCEDURE:
%
% NOTES:
% 	Implementation of the extended Hoshen-Kopelman algorithm to find clusters 
% 	from a set of points in N-dimmensional space and 1-dimensional time
%
%
% MODIFICATION HISTORY:
% 11/01/2011: Created by S.H.Lee, UC Berkeley
% 
%  Copyright (c) 2011 Sang-Hyuk Lee 
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

function [nCluster, idCluster] = fastsptcluster(posArr, t1Arr, t2Arr, rMax, tMax);

sz = size(posArr);
nArr = sz(1);

if rMax < 0 
   error('Maximum spatial distance should be a positive number');
end
if tMax < 0 
   error('Maximum temporal distance should be a positive number');
end

% Initialize array of neighbors
nNodeNextMax = 100;	% Guess enough number of neighbors  
nodeNext = zeros(nArr, nNodeNextMax);

nNodeNext = 0;

% Initialize node state. Every points are occupied.
nodeS = ones(nArr, 1);

% Sort in t1Arr
[t1Arr, idxT1] = sort(t1Arr);
t2Arr = t2Arr(idxT1);
posArr = posArr(idxT1,:);

% Find the first and the last indices of clusters with same value of t1Arr
firstIdxT1_cl = find(t1Arr ~= circshift(t1Arr, 1));
firstIdxT1_cl(1) = 1;        % First element always belongs to the first cluster
lastIdxT1_cl = circshift(firstIdxT1_cl, -1) - 1;
lastIdxT1_cl(end) = nArr; % Last element always belongs to the last cluster 

% Preallocate memory for arrays
firstIdxPreScreened =  zeros(nArr, 1);
lastIdxPreScreened =  zeros(nArr, 1);

nCl = numel(firstIdxT1_cl);

t1BoxHw = min(ceil(tMax), nCl);

id1_cl = [ones(t1BoxHw, 1);(1:nCl - t1BoxHw)'];
id2_cl = [(1+t1BoxHw:nCl)';nCl*ones(t1BoxHw,1)];

for i=1:nCl
   firstIdxPreScreened(firstIdxT1_cl(i):lastIdxT1_cl(i)) = firstIdxT1_cl(id1_cl(i));
   lastIdxPreScreened(firstIdxT1_cl(i):lastIdxT1_cl(i)) = lastIdxT1_cl(id2_cl(i));
end

for i=1:nArr
   % First, time adjacency 
   w = firstIdxPreScreened(i):lastIdxPreScreened(i);
   ww  =  find( t2Arr(i) - t1Arr(w) <= tMax ... 
        & t2Arr(w) - t1Arr(i) <= tMax );  
   w=w(ww);
   % Second, space adjacency
   currPosArr = posArr(i,:);
   ww = sum((posArr(w,:) - currPosArr(ones(numel(w), 1),:)).^2, 2) <= rMax^2;

   w = w(ww);

   w = w(w~=i); % Exclude itself from the list of its neighbors
   if (~isempty(w))
      nodeNext(i,1:numel(w))=w;
      nNodeNext = max(nNodeNext, numel(w));
   end
end

%fprintf('nNodeNext: %d \n',nNodeNext)

if nNodeNext < nNodeNextMax
   nodeNext = nodeNext(:,1:nNodeNext);
end

% Extended Hoshen-Kopelman clustering
[nCluster, idCluster] = fasthknodes(nodeS, nodeNext, 1);

% Reverse indexing to the unsorted array
idCluster(idxT1) = idCluster;

