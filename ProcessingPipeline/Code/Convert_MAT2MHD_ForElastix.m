%% inport image files and convert to .mhd
function [] = Convert_MAT2MHD_ForElastix(image, outpath, varargin)

% load([matpath filename]);
for ii=1:length(varargin)
    
    if isequal(varargin{ii},'outputname') 
        outputname=varargin{ii+1};
    end
    
    if isequal(varargin{ii},'convert_factor') 
        convert_factor=varargin{ii+1};
    end
    
    if isequal(varargin{ii},'sp') 
        sp=varargin{ii+1};
    end
    if isequal(varargin{ii},'orig')
        orig=varargin{ii+1};
    end
    if isequal(varargin{ii},'orient')
        orient=varargin{ii+1};
    end
    
end
%}
if ~exist('outputname','var')
    outputname = 'Output';
end
if ~exist('convert_factor','var')
    convert_factor = 1;
end
if ~exist('sp','var')
    sp = [1 1 1]';
end
if ~exist('orig','var')
    orig = [0 0 0]';
end
if ~exist('orient','var')
    orient = eye(3);
end


img = uint16(image);
mri = ImageType(size(img)', orig, sp, orient);
mri.data = img;

write_mhd([outpath outputname '.mhd'],mri);

% save([outpath outputname '_convertfactor.mat'], 'convert_factor')


