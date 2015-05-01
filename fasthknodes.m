%+
% NAME:
% 	fasthknodes
%
% PURPOSE:
% 	Accelerated extended Hoshen--Kopelmann clustering algorithm
%
% CATEGORY:
% 	Clustering
%
% CALLING SEQUENCE:
% 	[NumberOfClusters, NodeL]=fasthknodes(NodeS, NodeNext, OFlag)
%
% INPUTS:
%	NodeS: Occupancy state of nodes 
%	NodeNext: Neighboring nodes connected to each node 
%	OFlag: Flag for occupied nodes and links 
%
% OUTPUTS:
%	NumberOfClusters: Number of occupied clusters 
%	NodeL: Cluster labels of nodes 
%
% NOTES:
% 	Adapted from the original 'hknodes' code written by 'Ahmed AL-Futaisi and Tadeusz Patzek'
%
% % REFERENCE:
% 1. A. Al-Futaisi and T. W.. Patzek, 
%    "Extension of Hoshen-Kopelman algorithm to non-lattice environments," 
%     Physica A 321, 665-678 (2003).
% 
%
% MODIFICATION HISTORY:
% 11/07/2011: Created by S.H. Lee, University of California, Berkeley
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

function [NumberOfClusters, NodeL]=fasthknodes(NodeS, NodeNext, OFlag)
%===================================================== 
% 
% Adaptation of Hoshen--Kopelman cluster labeling algorithm for disordered 
% networks (nodes and links are placed at random in space)
% 
% Input arguments: 
% 
%	NodeS 		=Occupancy state of nodes 
%	NodeNext 	= Neighboring nodes connected to each node 
%	OFlag 		= Flag for occupied nodes and links (1 in this paper)
% 
% Output arguments: 
% 
%	NumberOfClusters= Number of occupied clusters 
%	NodeL 		= Cluster labels of nodes 
% 
% By: Ahmed AL-Futaisi and Tadeusz Patzek 
% 
%=========================================
% 
% STEP 1: READ THE DATA AND INITIALIZE THE OUTPUT 
% 
NumberOfNodes	= length(NodeS);	% Number of nodes in current network 
% 
% STEP 2: INITIALIZE THE HK ALGORITHM VARIABLES NodeL 
% 
NodeL=zeros(NumberOfNodes,1);	% Array to store cluster labels of nodes 
% 
% STEP 3: CREATE EMPTY ARRAY NodeLP AND START CLUSTER COUNTER 
% 
NodeLP = zeros(NumberOfNodes,1); % Array used for relabeling steps 
Cluster=0;	% Cluster counter
% 
% STEP 4: SCAN THE NETWORK NODES 
% 
for i=1:NumberOfNodes
% 
% Check if the node (Case 4c): 
%	1. has OFlag occupancy 
%	2. has NodeNext that have OFlag occupancy 
	N=find((NodeS(i)==OFlag).*(NodeS(nonzeros(NodeNext(i,:)))==OFlag));
       	if (~isempty(N))
	   %
	   % Define the occupancy status of NodeNext 
	   % 
	   Nodes=NodeNext(i,N); 
	   NodeNextL=NodeL(Nodes);
	   % 
	   % Case 4c i: No labeled neighbour 
	   % 
	   if any(NodeNextL)==0	% Start a new cluster
	      Cluster=Cluster+1; 
	      NodeL(i)=Cluster; 
	      NodeLP(Cluster)=Cluster;
	      % 
	      %  Case 4c ii: There exists a labeled neighbor 
	      % 
	   else % Put in the minimum labeling
	      % ====== Original ==================== 
	      %N=(NodeLP(nonzeros(NodeNextL))); 
	      %NodeLPmin=min(NodeLP(N)); 
	      %NodeL(i)=NodeLPmin; 
	      %NodeLP(N)=NodeLPmin; 
	      % ====================================
	      % 
	      % ===  Modification by Metzger et al Physica A 363 (2006) 558-560 ===
	      N = nonzeros(NodeNextL)';
	      for k = 1:length(N)
		 M = NodeLP(N(k));
		 while M < N(k)
		    N(k) = M;
		    M = NodeLP(N(k));
		 end
	      end
	      NodeLPmin = min(N);
	      NodeL(i) = NodeLPmin;
	      NodeLP(N) = NodeLPmin;
	      % ====================================
	   end
	% 
	%  This node is type 4b: 
	elseif NodeS(i)==OFlag
	   Cluster=Cluster+1; % Start a new cluster 
	   NodeL(i)=Cluster; 
	   NodeLP(Cluster)=Cluster; end
       	% 
	%  Skip nodes that are type 4a
end
% Resize NodeLP
NodeLP = NodeLP(1:Cluster);
%
% STEP 5A: CORRECT LABELS IN NodeLP RECURSIVELY
%
for i=1:length(NodeLP) 
   N=i; 
   while (NodeLP(N)<N)
      N=NodeLP(N); 
   end
   NodeLP(i)=N;
end
%
%	STEP 5B: RENUMBER LABELS IN NodeLP TO RUN SEQUENTIALLY
%

% 11/06/2011: Modified by S.H.Lee to speed up the calculation 
[uniqVal, id1, id2, idRev] = sortunique(NodeLP);
nUniq = numel(uniqVal);
newUniqNodeLP  =  (1:nUniq)';
w=find(uniqVal(:) ~= newUniqNodeLP);
for i=1:numel(w)
   NodeLP( idRev(id1(w(i)):id2(w(i))) ) = newUniqNodeLP(w(i));
end

% ================= Original  ================================
%NodeLP1=sort(NodeLP); 
%RelabL=NodeLP1(2:end).*(NodeLP1(2:end) > NodeLP1(1:end-1)); 
%RelabL=[NodeLP1(1), nonzeros(RelabL)'];
%for i=1:length(RelabL) 
%   NodeLP(find(NodeLP==RelabL(i)))=i;
%end
% ============================================================ 

% 
%	STEP 6: APPLY THE CORRECT LABELS TO THE ARRAYS NodeL
%


% 11/06/2011: Modified by S.H.Lee to speed up the calculation 
[uniqVal, id1, id2, idRev] =sortunique(NodeL);
if uniqVal(1) == 0 
   w = 2:numel(uniqVal);
else
   w = 1:numel(uniqVal);
end

if numel(w) ~= numel(NodeLP) 
   error('Error in the number of NodeLP');
end;

for i=1:numel(w)
   NodeL( idRev(id1(w(i)):id2(w(i))) ) =  NodeLP(i);
end


% ================= Original  ================================
%for i=1:length(NodeLP) 
%   NodeL(find(NodeL==i))=NodeLP(i);
%end
% ============================================================ 

% 
% 	RECORD NUMBER OF CLUSTERS
% 
NumberOfClusters=max(NodeL); 
%
return
