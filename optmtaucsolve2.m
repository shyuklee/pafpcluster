%+
% NAME:
% 	optmtaucsolve2
%
% PURPOSE:
% 	Approximate solution of optimal tc that achieves unbiased molecular counting by  
% 	balancing the over- and under-counting of blinking PA-FPs with double exponential dwell in 
% 	the dark states.  
%
% CATEGORY:
%
% CALLING SEQUENCE:
%	tc = optmtaucsolve2(kb, kd, kr1, kr2, alpha, tf, n)
%
% INPUTS:
% 	kb: photobleaching rate [sec]
% 	kd: dark states entry rate [sec]
% 	kr1: slow dark state exit rate [sec]
% 	kr2: fast dark state exit rate [sec] (kr2 > kr1)
% 	alpha: ratio kr2-population to kr1-population 
% 		p(t) = (kr1 * exp(-kr1*t) + alpha * kr2 * exp(-kr2*t)) / (1 + alpha)
% 	tf: The time within which all the PA-FPs get photoactivated with uniform photoactivation rate.
% 	    Fermi time when Fermi photoactivation scheme is employed.
% 	n: number of PA-FP molecules
%
% OUTPUTS:
% 	tc: optimal tc value
%
% REFERENCE:
% 1. S. Lee, J. Y. Shin, A. Lee and C. Bustamante, 
%    "Counting single photoactivatable fluorescent molecules by 
%     photoactivated localization microscopy (PALM),"
%     Proc. Natl. Acad. Sci. 109, 17436-17441 (2012).
%
%
% MODIFICATION HISTORY:
% 03/12/2012: Sang-Hyuk Lee, UC Berkeley.
% 12/27/2013: SHL. Documented. 
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

function tc = optmtaucsolve2(kb, kd, kr1, kr2, alpha, tf, n)

fun = @(x) (1-(1-x/tf)^n)/((exp(-kr1*x)+alpha*exp(-kr2*x))/(1+alpha)) - kd/kb*n/(n-1);

tc = fzero(fun, 0);
