%Inclass 14

%Work with the image stemcells_dapi.tif in this folder

% (1) Make a binary mask by thresholding as best you can
img = imread('stemcells_dapi.tif');
img_bw = img>300;
imshow(img_bw);

% (2) Try to separate touching objects using watershed. Use two different
% ways to define the basins. (A) With erosion of the mask (B) with a
% distance transform. Which works better in this case?

%% erosion
cc = bwconncomp(img_bw);
stats = regionprops(cc, 'Area');
area = [stats.Area];
fusedCandidates = area > mean(area) + std(area);
sublist = cc.PixelIdxList(fusedCandidates);
sublist = cat(1, sublist{:});
fusedMask = false(size(img_bw));
fusedMask(sublist) = 1;
s = round(sqrt(mean(area))./pi);
nucmin = imerode(fusedMask, strel('disk', s));
outside  = -imdilate(fusedMask, strel('disk', 100));
basin = imcomplement(bwdist(outside));
basin = imimposemin(basin, nucmin | outside);
pcolor(basin), shading flat;
L = watershed(basin);
newNuclearMask = L>1 | (img_bw - fusedMask);
imshow(newNuclearMask, 'InitialMagnification', 'fit')

%% distance transform
D = bwdist(~img_bw);
D = -D;
D(~img_bw) = -inf;
L = watershed(D);
L(~img_bw) = 0;
newNuclearMask1 = L>1;
imshow(newNuclearMask1, 'InitialMagnification', 'fit')
