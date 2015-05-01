%+
% NAME: 
% 	sortunique
%
% PURPOSE:
%  	Obtain indices of elements with an equal value of a sorted array 
%
%
% CALLING SEQUENCE:
%	[valUnique, idFirstSort, idLastSort, idReverse] = sortunique(A) 
%
%
% INPUTS:
% 	A: array of arbitrary dimension
% 	
%
% OUTPUTS:
% 	valUnique: unique sorted value
% 	idFirstSort(i): first index of the sorted array with the value equal to valUnique(i)  
% 	idLastSort(i): last index of the sorted array with the value equal to valUnique(i)  
% 	idReverse: A(idReverse) = sort(A)
%
%
% MODIFICATION HISTORY:
% 11/07/2011: Created by S.H.Lee, UC Berkeley
% 29/12/2013: SHL, Documented 

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


function [valUnique, idFirstSort, idLastSort, idReverse] = sortunique(A) 

[A idReverse] = sort(A(:)); 

idFirstSort = find(A ~= circshift(A, 1));
idFirstSort(1) = 1;        % First element always belongs to the first cluster
idLastSort = circshift(idFirstSort, -1) - 1;
idLastSort(end) = numel(A); % Last element always belongs to the last cluster 

valUnique = A(idFirstSort);


