function [scan1,scan2,scan3] = SizingTumors(scan1,scan2,scan3,voxVolume,imagedims)
%``````````````````````````````````````````````````````````````````````````
%Cellularity and volume
Tot1 = 0;
Tot2 = 0;
Tot3 = 0;

Tot1 = [sum(scan1.NTC,'all'),length(find(scan1.NTC>=1))];
Tot2 = [sum(scan2.NTC,'all'),length(find(scan2.NTC>=1))];
Tot3 = [sum(scan3.NTC,'all'),length(find(scan3.NTC>=1))];

scan1.TotCells = Tot1(1);
scan2.TotCells = Tot2(1);
scan3.TotCells = Tot3(1);

scan1.Volmm3 = Tot1(2)*voxVolume;
scan2.Volmm3 = Tot2(2)*voxVolume;
scan3.Volmm3 = Tot3(2)*voxVolume;

%``````````````````````````````````````````````````````````````````````````
%Longest axis
%{
%Using function Regionprops3 which is now a part of the MATLAB 2018 package
%The function calculates the longest axis of pixels in 3D
%As our voxels are not 1x1x1 need to regrid the data to get an accurate
%approximation of the longest axis using regionprops3

    dcesize = size(scan1.roi);
    [X3D,Y3D,Z3D] = meshgrid(1:dcesize(2),1:dcesize(1),1:dcesize(3));
    %Defining new 1x1x1mm grid
    NewGridSize = dcesize(1)*imagedims(1);
    NewGridSizeZ = dcesize(3)*imagedims(3);
    dx = (dcesize(2)-1)/NewGridSize;
    dy = (dcesize(1)-1)/NewGridSize;
    dz = 1;
    gridx = 1:dx:dcesize(2);
    gridy = 1:dy:dcesize(1);
    gridz = 1:dz:NewGridSizeZ;
    [xq,yq,zq] = meshgrid(gridx,gridy,gridz);

for bb = 1:3
    clear Z Znewgrid Znewgridmask;
    %Pulling out tumor ROI
    eval(['Z = scan' num2str(bb) '.roi;']); 
    %Creating a 3D grid of the original sizes

    %Initializing value to zero
    L.MajorAxisLength = 0; 

    %Only calculates longest axis if there is an ROI
    if(sum(sum(sum(Z)))>0)
        %Convert original ROI to the new 1x1x1mm grid
        Znewgrid = griddata(X3D,Y3D,Z3D,Z,xq,yq,zq,'linear');
        Znewgridmask = Znewgrid;
        Znewgridmask(Znewgridmask>0.01) = 1;
        Znewgridmask(isnan(Znewgridmask)==1) = 0;        
        L = RegionProps3(Znewgridmask,'MajorAxisLength');
    end
    eval(['scan' num2str(bb) '.LongestAxismm = ' num2str(L.MajorAxisLength) ';']);
end
%}





