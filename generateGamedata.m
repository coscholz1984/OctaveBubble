% BubbleTiles = imread('Bubbles.png');
Tilemap = imread('Tileset.png');%horzcat(cImg,BubbleTiles);

% Example of Bubble image to determine alpha mask for transparency
bImg = Tilemap(1:16,17:32,:);
alphaMask = (bImg(:,:,1) == 255) & (bImg(:,:,3) == 255);

% Convert Tilemap to index and get palette
[Tilemap, palette] = rgb2ind(Tilemap);

% Store Tiles in image array
bImages = [];
for iTile = 0:Nbubbles
  bImages(:,:,iTile+1) = Tilemap(:,(1:16)+iTile*16);
end

% make sure tiles are stored as uint8 (number of colors < 8bit)
bImages = uint8(bImages); 

save('gamedata.mat','bImages','palette','alphaMask');
