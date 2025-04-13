%% inport image files and convert to .mhd
function [] = Convert_MAT2MHD_SlicerFriendly(image, outpath, varargin)

for ii = 1:length(varargin)
    if isequal(varargin{ii},'outputname') 
        outputname=varargin{ii+1};
    end
    if isequal(varargin{ii},'convert_factor') 
        convert_factor=varargin{ii+1};
    end
end

if ~exist('outputname','var')
    outputname = 'Output';
end
if ~exist('convert_factor','var')
    convert_factor = 1;
end

img = uint16(image); 

%Introduce origin for the voxel space (Default is [0 0 0]')
orig = [0 0 0]';

%Introduce voxel spacing (Default for BrCa is [4/3 4/3 5]' mm)
sp = [4/3 4/3 5]';

%Introduce trasformation matrix (Default is [0 1 0; 1 0 0; 0 0 1])
% orient = eye(3);
orient = [0 1 0; 1 0 0; 0 0 1];

%Build mhd files
mri = ImageType(size(img)', orig, sp, orient);
mri.data = img;
write_mhd([outpath outputname '.mhd'],mri);

