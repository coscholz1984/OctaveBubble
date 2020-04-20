## Copyright (C) 2020 Christian
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} getConnectivity (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Christian <cnceo@CC-TP-T460>
## Created: 2020-04-19

function [DTtriangles, DTedges] = getConnectivity(NodeLocation)

% determine the connectivity of particles from the Delaunay triangulation
DTtriangles = delaunay(NodeLocation);
% remove all delaunay edges at the left and right border
for iDT = 1:size(DTtriangles,1)
  if (sum(NodeLocation(DTtriangles(iDT,:),2) == 8.5) > 1) || (sum(NodeLocation(DTtriangles(iDT,:),2) == 120.5) > 1)
    DTtriangles(iDT,:) = deal(NaN);
  end
end
DTtriangles(isnan(DTtriangles(:,1)),:) = [];
DTedges = [DTtriangles(:,[1 2]); DTtriangles(:,[2 3]); DTtriangles(:,[1 3])];
DTedges = sort(unique(DTedges,'rows'),2); % delete all duplicates and sort

end
