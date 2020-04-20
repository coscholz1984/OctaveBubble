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
## @deftypefn {} {@var{retval} =} findConnectedNodes (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Christian <cnceo@CC-TP-T460>
## Created: 2020-04-19

function connectedNodes = findConnectedNodes (DTedges, iNewParticles)
  connectedEdges = ismember(DTedges,iNewParticles);           % find all connections to the newly add particle
  connectedEdges = connectedEdges(:,1) | connectedEdges(:,2); % list all delaunay edges that include new particle
  connectedNodes = unique(DTedges(connectedEdges,:));         % list all unique nodes that connect to newly add particle
end
