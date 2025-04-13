%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Code simulates model for indiviual patients

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Intializing the workspace
clear all
close all

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Define path for computations
Directory = pwd;

addpath(Directory);
cd(Directory)

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Enter specific patient identifiers
patientID = 'testpatient';

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Designate start scan for simulation
startscan = 2;
%Designate comparison/target scan
endscan = 3;

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Read in patient information for image dimensions and scan times
fid = fopen([patientID '.txt']);
tline = fgetl(fid);
counter = 1;
while ischar(tline)
    if(counter == 1)
        imagedims = tline;
    elseif(counter == 2)
        times = tline;
    elseif(counter == 3)
        schedule = tline;
    elseif(counter == 4)
        slicebounds = tline;  
    elseif(counter == 5)
        GridSize = tline;    
    end    
    counter = counter + 1;    
    tline = fgetl(fid);
end
fclose(fid);
imagedims = str2num(imagedims);
times = str2num(times);
%Final time in days from intial scan to comparison scan
schedule = strsplit(schedule);
schedule = char(schedule);
scans = find(schedule=='S');
tf = sum(times(scans(startscan)+1:scans(endscan))); 
slicebounds = str2num(slicebounds);
%Specific slice top and bottom for model simulation efficiency
slicestart = slicebounds(1);
sliceend = slicebounds(2);
whichslices = slicestart:sliceend;
pslices = numel(whichslices);
n = str2num(GridSize);

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Load scan day info
load(['NativeX_' patientID '.mat']);
load(['NativeY_' patientID '.mat']);
load(['BreastMask_' patientID '.mat']);
if startscan == 1
    load(['NTC1_' patientID '.mat']);
    load(['Tissues1_' patientID '.mat']);
    z = NTC1(:,:,whichslices);
    tissues = Tissues1(:,:,whichslices);
elseif startscan == 2
    load(['NTC2_' patientID '.mat']);  
    load(['Tissues2_' patientID '.mat']);
    z = NTC2(:,:,whichslices);
    tissues = Tissues2(:,:,whichslices);
elseif startscan == 3
    load(['NTC3_' patientID '.mat']);    
    load(['Tissues3_' patientID '.mat']);
    z = NTC3(:,:,whichslices);
    tissues = Tissues3(:,:,whichslices);
else
    load(['NTC4_' patientID '.mat']); 
    load(['Tissues4_' patientID '.mat']);
    z = NTC4(:,:,whichslices);
    tissues = Tissues4(:,:,whichslices);
end
bw = BreastMask(:,:,whichslices);

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Coarsening grid for computational efficiency
dx = (max(max(X))-min(min(X)))/(str2num(GridSize)/2);
dy = (max(max(Y))-min(min(Y)))/(str2num(GridSize)/2);
[xq,yq] = meshgrid(min(min(X)):dx:max(max(X)),min(min(Y)):dy:max(max(Y)));
gridsize = size(xq);
n = gridsize(1);

%Interpolate points using cubic interpolator for each slice
for kk = 1:size(bw,3)
    Z(:,:,kk) = griddata(X,Y,z(:,:,kk),xq,yq,'cubic');
    Tissues(:,:,kk) = griddata(X,Y,tissues(:,:,kk),xq,yq,'cubic');
    BW(:,:,kk) = griddata(X,Y,bw(:,:,kk),xq,yq,'cubic');
end
Z(Z<0) = 0;
Tissues(Tissues<0) = 0;
BW(BW<0) = 0;
BW = round(BW);
Tissues = round(Tissues);

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Setting time vectors and identifying boundaries
dt = 0.25;
%Make the new time vector for the simulations
tvec = 0:dt:tf;
timesteps = size(tvec,2);

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Creating boundaries around the breast mask and also for the biopsy marker
%Boundary about the breast
[BCF] = Boundaries3DUT(n,BW);

dx = imagedims(1);
dy = imagedims(2);
dz = imagedims(3);

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Carrying capacity calculations for max possible per voxel size
TCvolume = 4189e-9;
packingdensity = 0.7405;
voxelvolume = prod(imagedims);
tHeta = packingdensity*voxelvolume/TCvolume;

CarryCaps = ones(n,n,pslices);
CarryCaps = CarryCaps*tHeta;
CarryCaps(CarryCaps<1) = 1;
countThreshold = 0.25*max(max(max(CarryCaps)));

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Defining values for mechanical coupling of tissue properties to diffusion
%Shear modulus and Hooke's law matrix
nU = 0.45; 

%Convert elements for the tumor tissue
E = 2e3*ones(size(Tissues));            %Adipose
E(Tissues==2) = 2*2e3;                  %Fibroglandular
ScaledZ = Z/tHeta;
E(ScaledZ>0.1) = 10*2e3;                %Tumor
clear ScaledZ 

%Calculate shear modulus 
G = E/2*(1-nU);
clear E

%Derivatives of G
[DiffyxG,DiffyyG,DiffyzG] = Diffy3D(dx,dy,dz,BCF,G);
DiffyG = zeros(n,n,pslices,3);
DiffyG(:,:,:,1) = DiffyxG;
DiffyG(:,:,:,2) = DiffyyG;
DiffyG(:,:,:,3) = DiffyzG;
clear DiffyxG DiffyyG DiffyzG

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Loading calibrated parameters
DMatrix = ones(n,n,pslices);

params = load(['params_' patientID '.txt']);
DMatrix = params(1)*DMatrix;
kMatrix = reshape(params(2:end),[n,n,pslices]);

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Run the model
[TCcoarse] = NtcFDmech3DUTLogisticOnly(Z,imagedims,DMatrix,dt,kMatrix,...
                                      CarryCaps,BCF,timesteps,G,DiffyG,nU);    


%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Interpolate results back to original grid
for kk = 1:size(bw,3)
    TC(:,:,kk) = griddata(xq,yq,TCcoarse(:,:,kk),X,Y,'cubic');
end
TCthresh = TC;
TCthresh(isnan(TCthresh)) = 0;
TCthresh(TCthresh<countThreshold) = 0;
n = size(TCthresh,1);

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Calculating total cellularity and volume
TotCellsTh = sum(TCthresh,'all');
NodesTh = length(find(TCthresh>=1));

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Calculating the longest axis for the result
%Creating a 3D grid of the original sizes
matrixsize = size(TCthresh);
%Creating a mask
Z = TCthresh;
Z(Z>0) = 1;
Z(Z<0) = 0;
[X3D,Y3D,Z3D] = meshgrid(1:matrixsize(2),1:matrixsize(1),1:matrixsize(3));
%Defining new 1x1x1mm grid
NewGridSize = matrixsize(1)*imagedims(1);
NewGridSizeZ = matrixsize(3)*imagedims(3);
dx = (matrixsize(2)-1)/NewGridSize;
dy = (matrixsize(1)-1)/NewGridSize;
dz = 1;
gridx = 1:dx:matrixsize(2);
gridy = 1:dy:matrixsize(1);
gridz = 1:dz:NewGridSizeZ;
[xq,yq,zq] = meshgrid(gridx,gridy,gridz);

%Initializing value to zero
L.MajorAxisLength = 0; 

%Only calculates longest axis if there is an ROI
if(sum(sum(sum(Z)))>0)
    %Convert original ROI to the new 1x1x1mm grid
    Znewgrid = griddata(X3D,Y3D,Z3D,Z,xq,yq,zq,'linear');
    Znewgridmask = Znewgrid;
    Znewgridmask(Znewgridmask>0.01) = 1;
    Znewgridmask(isnan(Znewgridmask)) = 0;        
    L = regionprops3(Znewgridmask,'MajorAxisLength');
end    

cd(Directory)
beep

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end of file