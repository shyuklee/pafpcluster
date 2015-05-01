%+
% NAME:
% 	optmtaucsolve2_dendra
%
% PURPOSE:
% 	Approximate solution of optimal tc using in vitro rates of Dendra2 
%
% CATEGORY:
%
% CALLING SEQUENCE:
%	tc = optmtaucsolve2_dendra(nMol, tF)
%
% INPUTS:
% 	nMol: (M) array of number of molecules 
% 	tF: (N) array of Fermi time [min] 
%
% OUTPUTS:
% 	tc: (M, N) optimal tc value [sec]
%
% REFERENCE:
% 1. S. Lee, J. Y. Shin, A. Lee and C. Bustamante, 
%    "Counting single photoactivatable fluorescent molecules by 
%     photoactivated localization microscopy (PALM),"
%     Proc. Natl. Acad. Sci. 109, 17436-17441 (2012).
%
% MODIFICATION HISTORY:
% 03/12/2012: Sang-Hyuk Lee, UC Berkeley.
% 12/27/2013: SHL. Documented. 
%
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

function [tc, nMol, tF] = optmtaucsolve2_dendra(nMol, tF)

% Convert nMol to a column vector
nMol = nMol(:);
M = numel(nMol);

% Convert tF to a row vector
tF = tF(:)';
N = numel(tF);

% Make a grid
nMol = nMol(:,ones(1, N));
tF = tF(ones(M, 1), :);

% Dendra2 rates in vitro
kb = 16.6;
kd = 3.2;
kr1 = 1.6;
kr2 = 18.0;
alpha = 3.2;

% Approximate solver
for n = 1:N
   for m = 1:M
      tc(m, n) = optmtaucsolve2(kb, kd, kr1, kr2, alpha, 60*tF(m,n), nMol(m,n));
   end
end

