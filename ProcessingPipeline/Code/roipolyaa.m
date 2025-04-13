function mask=roipolyaa(xmin, xmax, ymin, ymax, XV, YV, sd)
%ROIPOLYAA 1 for all points inside polygon, zero for points outside, and
%   varying degrees between 1 and 0 for points on the boundary.
%   Usage is as for ROIPOLY in the imaging toolbox, with addition of
%   anti-aliasing as described at:
%   http://blogs.mathworks.com/steve/2007/02/06/antialiased-polygon-scan-conversion/
%
%   mask=ROIPOLYAA(xmin, xmax, ymin, ymax, XV, YV, sd) returns a matrix
%   MASK covering the region from (XMIN, YMIN) to (XMAX, YMAX).
%   MASK(p, q) = 1 if the point (p, q) is fully inside the polygon whose
%   vertices are specified by the vectors XV and YV. MASK(p, q) = 0 if the
%   point is fully outside the polygon. Points which are on the boundary
%   are weighted from 0 to 1 using the method discussed in the above URL,
%   and using a sub-pixel block dimension of SD.
%
%   SD must be an odd number (so that it can be centered on the element),
%   and larger than 1 (otherwise just use INPOLYGON directly). Higher
%   values for SD dramatically increase processing time.
%
%   All polygon types supported by INPOLYGON will work with ROIPOLYAA.
%
%   Example
%      XV = rand(4, 1) * 80 + 10;
%      YV = rand(4, 1) * 80 + 10;
%      mask = (1, 100, 1, 100, XV, YV, 5)
%
%   See also inpolygon, roipoly
%
%   Copyright 2014 Andrew Simpson


%% Sanity check the inputs
if sd <3 || mod(sd, 2) ~= 1
    error('Sub-pixel block dimension SD must be an even number over 1');
end
if xmax <= xmin ...
        || ymax <= ymin
    error('XMIN must be below XMAX, YMIN below YMAX');
end
if ~isequal(floor([xmin xmax ymin ymax]), [xmin xmax ymin ymax])
    error('XMIN, XMAX, YMIN, YMAX must be integers');
end

%% All data appears good, proceed

win=floor(sd/2); % window size: no of elements on each side of centre

% Create a mask of size [(xmax-xmin)*sd, (ymax-ymin)*sd] with all the
% vectors scaled up by sd. This is equivalent to sub-dividing each element
% into a block of size sd*sd
[xv, yv]=ndgrid(ymin*sd-win:ymax*sd+win, xmin*sd-win:xmax*sd+win);
maskL=double(inpolygon(xv, yv, YV*sd, XV*sd));

%% Create the mask by taking the mean of each block of sd*sd in maskL

% Pre-assign memory in mask variable (this makes no difference but is good
% practice anyway)
maskT=zeros(sd, numel(maskL)/sd);

% Split maskL down sd rows at a time, to allow 'mean' to be used better
for i=1:sd:size(maskL, 1)
    s=((i-1)/sd)*size(maskL, 2)+1;  % start column of temp mask
    e=((i-1)/sd+1)*size(maskL, 2);  % end column of temp mask
    maskT(:, s:e)=maskL(i:i+sd-1, :);  % build temp mask
end

% take mean of every column of temp mask.
maskT=mean(maskT);
% new temp mask is 1xn, where each sdxsd group from maskL is now 1xsd
% reshape temp mask so that it is sdxn, and each group from maskL is sdx1
maskT=reshape(maskT, sd, []);
% take the mean of each sdx1 column in temp mask
% m4 has all correct values for mask, just arranged in a single row
maskT=mean(maskT);
% reshape into the correct mask shape.
maskT=reshape(maskT, size(maskL,2)/sd, []);
mask=maskT';


%% Comment out the 'return' below to plot the results of the mask
return;

% If outlines is set to 1, print both the polygon outline and the point of
% interest
outlines=1;                                                    %#ok<*UNRCH>

% Point of interest for plots - one vertex of triangle
x=round(XV(3));
y=round(YV(3));

% For plotting, ensure coords are stored in first dimension
XVs=shiftdim(XV);
YVs=shiftdim(YV);

%% Plot the large mask created with inpolygon
figure; hold on;
image(((xmin*sd-2):(xmax*sd+2))/sd, ((ymin*sd-2):(ymax*sd+2))/sd, double(maskL)*100);
if outlines
    plot([XVs; XVs(1)], [YVs; YVs(1)], 'rx-');

    % Plot a square around the point of interest
    xs=(x*sd-sd/2)/sd; xe=(x*sd+sd/2)/sd;
    ys=(y*sd-sd/2)/sd; ye=(y*sd+sd/2)/sd;
    plot([xs xs xe xe xs], [ys ye ye ys ys], 'g-'); 
end

axis image;
colormap(gray);
title('Large mask');

%% Plot the finished mask
figure; hold on;
image(xmin:xmax, ymin:ymax, mask*100);

if outlines
    plot([XVs; XVs(1)], [YVs; YVs(1)], 'rx-');
    % Plot an X at the point of interest (the average value of the box in
    % the previous graph)
    plot(x, y, 'gx');
end

axis image;
colormap(gray);
title(['Mask, antialiased with block size ' num2str(sd) 'x' num2str(sd)]);
imwrite(mask, 'test1.png');

%% Finally, display what it would look like without anti-aliasing
[xv1, yv1]=ndgrid(ymin:ymax, xmin:xmax);
maskB=inpolygon(xv1, yv1, YV, XV);

figure; hold on;
image(xmin:xmax, ymin:ymax, maskB*100);

if outlines
    plot([XVs; XVs(1)], [YVs; YVs(1)], 'rx-');
    plot(x, y, 'gx');
end

axis image;
colormap(gray);
title('Mask, no anti-aliasing');
imwrite(maskB, 'test2.png');
