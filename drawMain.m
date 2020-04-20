% =================== Octave Puzzle Bubble  =================== %
% This is an implementation of the puzzle bubble game, also     %
% known from various clones, such as frozen bubble, but written %
% from scratch in GNU Octave. The game runs rather slow but     %
% kind of playable on modern hardware (~ May 2020)              %
%
% This code is being released for learning purposes and in the  %
% that it is usefull an funny. However WITHOUT ANY WARRANTY;    % 
% without even the implied warranty of MERCHANTABILITY or       % 
% FITNESS FOR A PARTICULAR PURPOSE.                             %
%
% The code is being licensed under the                          %
% GNU General Public license v3.0                               %
% Graphics and data used in this example are licensed under the %
% Creative Commons Attribution 4.0 International license        %
% (CC BY 4.0)                                                   %
%
% Author: Christian Scholz (2020)

% %%%%%%%%%%%%%%%%%%%%%%%%% Preamble %%%%%%%%%%%%%%%%%%%%%%%%% %
% | In this section we load the required packages and set    | %
% | some core variables                                      | %

% Load image package, which is required for color palette usage
pkg load image

% DBUG and ASSIST mode
DDBUG = false;
ASSIST = true;

% This loads the game graphics
load('gamedata.mat');

% %%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%% %
% | In this section we initialize some global variables and    | %
% | graphics for the game                                      | %

% Variables to store node location and type
NodeLocation = [];
NodeType = [];

% Height of playfield in numbers of bubbles
iTileN = 11;
% Number of bubble colors (must be <= 8)
Nbubbles = 4;
if Nbubbles > 8
  disp('Too many number of bubble colors.')
  return;
end

% angle for shooting the bubble
sAngle = 90;  

% true if a new particle has been placed
setNew = true;              
                
% generate background
Canvas = repmat(bImages(:,:,1),[iTileN 8]);
CanvasBG = Canvas;
% Plot background, initial particles and connections
if DDBUG
  figure;
end
% initialize and plot bubbles
for iTile = 1:iTileN
    if mod(iTile,2) == 1
      jTileN = 8;
      xshift = 0;
    else
      jTileN = 7;
      xshift = 8;
    end
    for jTile = 1:jTileN
        NodeLocation = [NodeLocation; (iTile-1)*14+8.5,(jTile-1)*16+8.5+xshift];
        if iTile < iTileN-1
          NodeType = [NodeType; randi(Nbubbles+1)];
        else
          NodeType = [NodeType; 1]; % last two rows should be empty
        end
        tmpImg = Canvas((iTile-1)*14+(1:16),(jTile-1)*16+(1:16)+xshift);
        bImgtmp = bImages(:,:,NodeType(end));%bImg(:,:,permutations(NodeType(end),:));
        if NodeType(end) > 1
          tmpImg(~alphaMask) = bImgtmp(~alphaMask);
        end
        Canvas((iTile-1)*14+(1:16),(jTile-1)*16+(1:16)+xshift) = tmpImg;
    end
end 

if DDBUG
  image(Canvas);
  colormap(palette);
  axis image;
end

% get connected nodes from modified delaunay triangulation
[DTtriangles, DTedges] = getConnectivity(NodeLocation);
% plot delaunay triangulation for visualization
if DDBUG
  hold on;
  triplot(DTtriangles,NodeLocation(:,2),NodeLocation(:,1)); 
  hold off;
end

% %%%%%%%%%%%%%%%%%%%%%%% Main Game Code %%%%%%%%%%%%%%%%%%%%%%% %
% | Here we start with the main code, draw the canvas and      | %
% | run the main game loop                                     | %

% draw the main game window
hFigure = figure;
haxes = gca;
colormap(palette);

% loop variable to hold the current NodeTypes during the game loop
newNodeType = NodeType;

isPlaced = true;    % set true if a new particle has been added
iNewParticles = 84; % dummy value for where a new particle is added

% main game loop
while true
  
  % chose type of new particle to add, make sure it is only of a type that is still 
  % existing on the field
  if setNew
    nTypes = unique(newNodeType);
    newParticleType = nTypes(randi(numel(nTypes)-1)+1);
  end
  
  % %%%%%% Find connected clusters of particles %%%%%% %
  
  % Cluster detection based on adding a new particle and search only connections to 
  % the newly added particle.
  if isPlaced
    connectedEqualNodes = [];
    while true
      connectedNodes = findConnectedNodes(DTedges, iNewParticles);
      if isempty(connectedNodes)
        if DDBUG 
          disp('No connectedNodes found');
        end
        break;
      end
      % keep only particles with equal type
      connectedEqualNodes = connectedNodes(ismember(newNodeType(connectedNodes),newNodeType(iNewParticles)));
      if numel(iNewParticles) == numel(connectedEqualNodes)
        break;  % if no new particles of equal type are found, exit the loop
      else
        iNewParticles = connectedEqualNodes;
      end
    end
  end
  
  % %%%%%% Delete particles if cluster >= 3 %%%%%% %
  
  % This will delete (set to type 1) the particles which are connected to the newly set 
  % node, if there is a total number of more than 2. Be warned that this will not work 
  % as intended if more than a single initial particle is specified, because it could happen 
  % that two or more clusters of identical type are found which could result in a total 
  % of more than 2, but they are actually not connected. In a proper game, this situation 
  % should never occur.
  if isPlaced
    % check if there are at least three non-empty equal connected particles
    if (numel(newNodeType(connectedEqualNodes)) > 2) && (range(newNodeType(connectedEqualNodes)) == 0) && (newNodeType(connectedEqualNodes(1)) > 1)
      if DDBUG
        disp('Cluster found! Removing it.')
      end
      newNodeType(connectedEqualNodes) = 1;
    end
  end
  
  % %%%%%% Remove floating bubbles and determine accesible region %%%%%% %
  
  % Here we search freely floating particles to remove them.
  % For this we can restrict to distinguish only empty nodes from occupied nodes.
  % A particle is floating if there is no path to any of the occupied nodes in the
  % first row. We loop through connections using a burning algorithm:
  % 1st pass: Find all occupied particles in top row
  % 2nd pass: Find all occpuied particles connected to occupied particles in top row
  % 3rd pass: Find all occupied particles connected to these etc....
  % repeat this until no more connections to occupied particles can be found.
  if isPlaced
    iUpperParticles = 1:8;
    iUpperParticles = iUpperParticles(newNodeType(iUpperParticles) ~= 1);
    while true
      connectedParticles = findConnectedNodes(DTedges, iUpperParticles);
      connectedEqualParticles = connectedParticles(newNodeType(connectedParticles) ~= 1); % keep only non-empty sites
      if numel(iUpperParticles) == numel(connectedEqualParticles)
        if DDBUG
          disp(['Hangig cluster has ', num2str(numel(iUpperParticles)), ' particles']);
        end
        break;
      else
        iUpperParticles = connectedEqualParticles;
      end
    end
    
    % delete all particles floating in air
    newNodeType(setdiff(1:83,iUpperParticles)) = 1;
    
    % search for lower accessible (empty) nodes
    iLowerParticles = 76:83;
    iLowerParticles = iLowerParticles(newNodeType(iLowerParticles) == 1);
    while true
      connectedParticles = findConnectedNodes(DTedges, iLowerParticles);
      connectedEqualParticles = connectedParticles(newNodeType(connectedParticles) == 1); % keep only empty sites
      if numel(iLowerParticles) == numel(connectedEqualParticles)
        if DDBUG
          disp(['Lower cluster has ', num2str(numel(iLowerParticles)), ' particles']);
        end
        break;
      else
        iLowerParticles = connectedEqualParticles;
      end
    end
    
    % search unoccupied ceiling particles
    iCeilParticles = 1:8;
    iCeilParticles = iCeilParticles(newNodeType(iCeilParticles) == 1);
    CeilLocations = NodeLocation(iCeilParticles,[2 1]);
    
    % find all nodes on the boundary between empty and occupied space
    % get list of occupied nodes
    occupiedNodes = find(newNodeType ~= 1);
    boundaryNodes = unique(DTedges(xor(ismember(DTedges(:,1),occupiedNodes),ismember(DTedges(:,2),occupiedNodes)),:));
    
    % get coordinates to check when the line first hits an occupied boundary particle 
    checkNodes = boundaryNodes(newNodeType(boundaryNodes) ~= 1);
    checkNodesLocation = NodeLocation(checkNodes,[2 1]); % x and y are fliped
    checkNodesEmpty = intersect(iLowerParticles,boundaryNodes(newNodeType(boundaryNodes) == 1));
    checkNodesEmptyLocation = NodeLocation(checkNodesEmpty,[2 1]); % x and y are fliped
  end
  
  % %%%%%% Find where to place new particle via ray tracing %%%%%% %
  
  % do the raytracing from the initial particle to the next occupied node
  % first, check if ray hits an occupied particles, i.e. nearest occupied particle 
  % to ray has a distance of less than a particle diameter (16 px)
  currentX = 64.5;
  currentY = 166.5;
  currentAngle = sAngle;
  
  lineX = currentX;
  lineY = currentY;
  lineAngle = currentAngle;
  bHit = false;
  
  dxPhit = NaN;
  dyPhit = NaN;
  
  % loop until either a bubble or the ceiling has been hit
  while (currentY > 0) && ~bHit
    lowlimit = 180-90-atand((128-currentX)/currentY);
    uplimit = 90+atand(currentX/currentY);
    if (currentAngle > lowlimit) && (currentAngle < uplimit)
      nextY = 0;
      nextX = currentX + currentY/tand(currentAngle);
    elseif currentAngle < lowlimit
      nextX = 128;
      nextY = currentY-(128-currentX)*tand(currentAngle);
    elseif currentAngle > uplimit
      nextX = 0;
      nextY = currentY+currentX*tand(currentAngle);
    end

    % distances between nodes and line
    odistances = ((nextY-currentY)*checkNodesLocation(:,1) - (nextX-currentX)*checkNodesLocation(:,2) + nextX*currentY - nextY*currentX)./sqrt((nextX-currentX)^2+(nextY-currentY)^2);
    ldistances = abs(odistances);
    % are there any distances below or equal bubble diameter 
    bdistances = (ldistances <= 16);
    % if line hits any particles check which is closest to currentX, currentY
    if any(bdistances)
      tmpNodesLocation = checkNodesLocation(bdistances,:);
      tmpdistances = ldistances(bdistances);
      tmpodistances = odistances(bdistances);
      tmpNodesIdx = checkNodes(bdistances);
      [~,didx] = min(sqrt((checkNodesLocation(bdistances,1)-currentX).^2 + (checkNodesLocation(bdistances,2)-currentY).^2));
      if DDBUG
        disp(['Found hit at ',num2str(tmpNodesLocation(didx,:))]);
      end
      bHit = true;
      nextX = tmpNodesLocation(didx,1);
      nextY = tmpNodesLocation(didx,2);
      hitParticle = tmpNodesIdx(didx);
      % now that we found the hit particle, we need to figure out where it is hit,  
      % and which empty node the new particle to place.
      % get the unit vector in the direction of propagation (pex, pey)
      [pex pey] = pol2cart(deg2rad(currentAngle),1);
      % get the -90 deg rotated vector to that (dex, dey)
      dex = pey;
      dey = -pex;
      % if we add -r * (pex,pey) and distance * (dex, dey) we should get the vector
      % between the centers of the colliding circles. Half of this is the point of contact.
      tmpDeltaLength = sqrt(16^2-tmpodistances(didx)^2);
      dxPhit = .5 * (-tmpDeltaLength * pex + tmpodistances(didx) * dex);
      dyPhit = .5 * (-tmpDeltaLength * pey + tmpodistances(didx) * dey);
      % calculate point of collision
      colx = tmpNodesLocation(didx,1) + dxPhit;
      coly = tmpNodesLocation(didx,2) - dyPhit;
      % find the nearest unoccipied boundary node to (colx,coly) and place particle there
      [~,uidx] = min(sqrt((checkNodesEmptyLocation(:,1) - colx).^2 +  (checkNodesEmptyLocation(:,2) - coly).^2));
      newParticleNode = checkNodesEmpty(uidx);
      % end of line, where particle hit another boundary particle
      nextX = tmpNodesLocation(didx,1) + 2*dxPhit;
      nextY = tmpNodesLocation(didx,2) - 2*dyPhit;
    elseif nextY <= 0
      disp('Hit ceiling');
      bHit = true;
      [minCeil,cidx] = min(sqrt((CeilLocations(:,1) - nextX).^2 +  (CeilLocations(:,2) - nextY).^2));
      [minEmpty,uidx] = min(sqrt((checkNodesEmptyLocation(:,1) - nextX).^2 +  (checkNodesEmptyLocation(:,2) - nextY).^2));
      if minEmpty < minCeil
        newParticleNode = checkNodesEmpty(uidx);
      else
        newParticleNode = iCeilParticles(cidx);
      end
    end
  
    currentX = nextX;
    currentY = nextY;
    currentAngle = 180 - currentAngle;
    lineX = [lineX; currentX];
    lineY = [lineY; currentY];
    lineAngle = [lineAngle; currentAngle];
  end
  
  % %%%%%% Draw the current frame %%%%%% %
  
  if isPlaced  
    pause(.2);
    Canvas = CanvasBG; 
    cNodes = 1;
    for iTile = 1:11
        if mod(iTile,2) == 1
          jTileN = 8;
          xshift = 0;
        else
          jTileN = 7;
          xshift = 8;
        end
        for jTile = 1:jTileN
            tmpImg = Canvas((iTile-1)*14+(1:16),(jTile-1)*16+(1:16)+xshift);
            bImgtmp = bImages(:,:,newNodeType(cNodes));
            if newNodeType(cNodes) > 1
              tmpImg(~alphaMask) = bImgtmp(~alphaMask);
            end
            Canvas((iTile-1)*14+(1:16),(jTile-1)*16+(1:16)+xshift) = tmpImg;
            cNodes = cNodes + 1;
        end
    end
    % plot new particle to add
    iTile = 12;
    jTile = 4;
    yshift = 4;
    xshift = 8;
    tmpImg = Canvas((iTile-1)*14+(1:16)+yshift,(jTile-1)*16+(1:16)+xshift);
    bImgtmp = bImages(:,:,newParticleType);
    tmpImg(~alphaMask) = bImgtmp(~alphaMask);
    Canvas((iTile-1)*14+(1:16)+yshift,(jTile-1)*16+(1:16)+xshift) = tmpImg;
  end

  image(haxes,Canvas);
  axis image;
  set(haxes,'YTickLabel',[]);
  set(haxes,'XTickLabel',[]);
  hold on;
  quiver(haxes,64.5,166.5,15*cosd(sAngle),-15*sind(sAngle),'-o','Color',[1.0 1.0 0.0],'LineWidth',3,'filled');
  if ASSIST
    quiver(haxes,tmpNodesLocation(didx,1),tmpNodesLocation(didx,2),tmpodistances(didx) * dex,-tmpodistances(didx) * dey,'-o','Color',[0.0 0.3 1.0],'LineWidth',3,'filled');
    plot(haxes,lineX, lineY, 'LineStyle', ':', 'LineWidth',5, 'Color', [.0 .4 .7]);
    plot(haxes,tmpNodesLocation(didx,1) + dxPhit, tmpNodesLocation(didx,2) - dyPhit,'o','MarkerSize',7);
    plot(haxes,NodeLocation(newParticleNode,2), NodeLocation(newParticleNode,1),'o','MarkerSize',7,'Color','r');
  end
  line(haxes,[0 129],[144 144], 'LineStyle', '-', 'LineWidth',3, 'Color', [1.0 .4 .2]);
  hold off;
  drawnow;
  
  % %%%%%% Determine if game has been won or lost %%%%%% %
  
  if sum(newNodeType(76:83) ~= 1) > 0
    disp('Game Over');
    hold on;
    text(128/2,146/2,"Game Over",'color',[1 1 1],'fontsize',32,'horizontalalignment',"center",'verticalalignment',"middle");
    hold off;
    break;
  end
  
  % end game if all particles are cleared
  if sum(newNodeType) == 83
    disp('Congratulations!');
    hold on;
    text(128/2,146/2,"Congratulations!",'color',[1 1 1],'fontsize',32,'horizontalalignment',"center",'verticalalignment',"middle");
    hold off;
    break;
  end
  
  % %%%%%% User Input %%%%%% %
  
  isPlaced = false;
  while ~isPlaced
    setNew = false;
    % record mouse click
    [x,y,buttons]=ginput(1);
    if DDBUG
      disp(['Mouse clicked on position ', num2str([x,y])]);
    end
    % Check if you clicked on a node
    % Calculate distance to nearest node
    if DDBUG
      disp(num2str(buttons));
    end
    if buttons == 97
      sAngle = sAngle + 3*0.8;
      if sAngle > 170
        sAngle = 170;
      end
      break;
    end
    if buttons == 100
      sAngle = sAngle - 3*0.8;
      if sAngle < 10
        sAngle = 10;
      end
      break;
    end
    if buttons == 32
      iNewParticles = newParticleNode;
      newNodeType(newParticleNode) = newParticleType;
      setNew = true;
      isPlaced = true;
      break;
    end
  end
  
  % %%%%%% Draw intermediate frame with all new particles %%%%%% %
  
  if isPlaced 
    Canvas = CanvasBG; 
    cNodes = 1;
    for iTile = 1:11
        if mod(iTile,2) == 1
          jTileN = 8;
          xshift = 0;
        else
          jTileN = 7;
          xshift = 8;
        end
        for jTile = 1:jTileN
            tmpImg = Canvas((iTile-1)*14+(1:16),(jTile-1)*16+(1:16)+xshift);
            bImgtmp = bImages(:,:,newNodeType(cNodes));
            if newNodeType(cNodes) > 1
              tmpImg(~alphaMask) = bImgtmp(~alphaMask);
            end
            Canvas((iTile-1)*14+(1:16),(jTile-1)*16+(1:16)+xshift) = tmpImg;
            cNodes = cNodes + 1;
        end
    end
    image(haxes,Canvas);
    axis image;
    hold on;
    line(haxes,[0 129],[144 144], 'LineStyle', '-', 'LineWidth',3, 'Color', [1.0 .4 .2]);
    hold off;
    drawnow;
  end

end
