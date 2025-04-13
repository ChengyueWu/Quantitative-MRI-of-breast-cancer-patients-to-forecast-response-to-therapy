function SetupRegFiles_ISPY2(subname, inpath, outpath, Visits,Visit_num_f, weight, varargin)
%All codes adapted from Chengyue Wu's function files for registration with
%Elastix using a Rigid intial registration followed by Bspline with a
%penalty for the tumor ROI to conserve size and shape of the tumor.

%Registration tends to perform better with blank slides
%"padding" each of the ends of the registration--otherwise parts of the
%image are lost (especially at the ends)

%! Registration cannot handle NaNs or Infs so they are removed from all
%maps
%``````````````````````````````````````````````````````````````````````````
% Load in data from pipeline

% inpath  = [Parentpath, subname,'/UnregisteredData/'];
% outpath = [Parentpath, subname,'/IntervisitRegistered/RegFiles/'];

if isfolder(outpath)==0
    mkdir(outpath)
end

%%%%
for ii=1:length(varargin)
    if isequal(varargin{ii},'YRange') 
        YRange=varargin{ii+1}; 
    end  
    if isequal(varargin{ii},'XRange') 
        XRange=varargin{ii+1};
    end 
    if isequal(varargin{ii},'ZRange') 
        ZRange=varargin{ii+1}; 
    end 
end
%}
%% Loading files that contain DCE data
DCE_multiV = cell(size(Visits));
ROI_multiV = cell(size(Visits));
ADC_multiV = cell(size(Visits));
for v = Visits
    filename = [subname 'T',num2str(v-1), '_DCEtc_3DRigidReg.mat']; 
    DCE_multiV{v} = load([inpath filename]);
    
    %%% CHANGE HERE %%%
%     filename = [subname 'T',num2str(v-1), '_DCEmask_FCM_3DRigidReg.mat']; 
%     ROI_multiV{v} = load([inpath filename]);
    filename = [subname 'T',num2str(v-1), '_ParaMap_DCEmask_align.mat']; 
    ROI_multiV{v} = load([inpath filename]);

    %%% CHANGE HERE %%%
%     filename = [subname 'T',num2str(v-1), '_ParaMap_ADC_3DRigidReg.mat']; 
%     ADC_multiV{v} = load([inpath filename]);
    filename = [subname 'T',num2str(v-1), '_Calculated_ADC_3DRigidReg.mat']; 
    ADC_multiV{v} = load([inpath filename]);
end

res_list = [DCE_multiV{1}.DCE_reg.resolution(1); ...
            DCE_multiV{2}.DCE_reg.resolution(1); ...
            DCE_multiV{3}.DCE_reg.resolution(1);...
            ];
scale_list = res_list./max(res_list);

%%%%% resample to make sure all images have the same resolution
% Downsample to visit with lowest resolution
    for v = Visits
        scale = scale_list(v);
        if scale ~= 1
            
            [Ysize_interp, Xsize_interp, Zsize_interp] = size(ROI_multiV{v}.ParaMap_DCEmask_align.Data);
            [X,Y,Z]    = meshgrid(1:Xsize_interp,1:Ysize_interp,1:Zsize_interp);
            [Xq,Yq,Zq] = meshgrid(1:1/scale:Xsize_interp, ...
                                  1:1/scale:Ysize_interp, ...
                                  1:Zsize_interp);         
            ROI_resample = interp3(X,Y,Z,double(ROI_multiV{v}.ParaMap_DCEmask_align.Data),Xq,Yq,Zq);
            %ADC_resample = interp3(X,Y,Z,double(ADC_multiV{v}.ParaMap_ADC_reg.Data),Xq,Yq,Zq);
            ADC_resample = interp3(X,Y,Z,double(ADC_multiV{v}.Calculated_ADC_reg),Xq,Yq,Zq);
            
            DCE_resample = zeros(size(Xq));
            for t = size(DCE_multiV{v}.DCE_reg.Data, 4)
                DCE_resample(:,:,:,t) = interp3(X,Y,Z,double(DCE_multiV{v}.DCE_reg.Data(:,:,:,t)),Xq,Yq,Zq);
            end
            
            %ADC_multiV{v}.ParaMap_ADC_reg.Data = ADC_resample;
            ADC_multiV{v}.Calculated_ADC_reg         = ADC_resample;
            ROI_multiV{v}.ParaMap_DCEmask_align.Data = ROI_resample;
            DCE_multiV{v}.DCE_reg.Data               = DCE_resample;
        end
    end
    
    
%%
if ~exist('YRange','var')
    Ysize  = [];
    for v = Visits
        Ysize = [Ysize; size(DCE_multiV{v}.DCE_reg.Data, 1)];
    end
    ysize_final = min(Ysize);
    YRange = [floor((Ysize(1)-ysize_final)/2) + (1:ysize_final); ...
              floor((Ysize(2)-ysize_final)/2) + (1:ysize_final); ...
              floor((Ysize(3)-ysize_final)/2) + (1:ysize_final)]; %repmat(1:min(Ysize), [length(Visits),1]);
end
if ~exist('XRange','var')
    Xsize = [];
    for v = Visits
        Xsize = [Xsize; size(DCE_multiV{v}.DCE_reg.Data, 2)];
    end
    xsize_final = min(Xsize);
    XRange = [floor((Xsize(1)-xsize_final)/2) + (1:xsize_final); ...
              floor((Xsize(2)-xsize_final)/2) + (1:xsize_final); ...
              floor((Xsize(3)-xsize_final)/2) + (1:xsize_final)]; %repmat(1:min(Xsize), [length(Visits),1]);
end
if ~exist('ZRange','var')
    Zsize = [];
    for v = Visits
        Zsize = [Zsize; size(DCE_multiV{v}.DCE_reg.Data, 3)];
    end
    zsize_final = min(Zsize);
    ZRange = [floor((Zsize(1)-zsize_final)/2) + (1:zsize_final); ...
              floor((Zsize(2)-zsize_final)/2) + (1:zsize_final); ...
              floor((Zsize(3)-zsize_final)/2) + (1:zsize_final)]; %repmat(1:min(Zsize), [length(Visits),1]);
end

%%
for v = Visits
    
    avgdce = mean(DCE_multiV{v}.DCE_reg.Data(YRange(v,:), XRange(v,:), ZRange(v,:),:),4);
    avgdce(isnan(avgdce)) = 0;
    avgdce(isinf(avgdce)) = 0;
    
    anatomical = zeros(size(avgdce));
    for mm = 1:size(anatomical,3)
        temp = avgdce(:,:,mm)/max(max(avgdce(:,:,mm)));
        temp(isnan(temp)) = 0;
        temp(isinf(temp)) = 0;
        anatomical(:,:,mm) = 350*adapthisteq(temp);
    end
    anatomical(isnan(anatomical)) = 0;
    anatomical(isinf(anatomical)) = 0;

    
    roi_orig = double(ROI_multiV{v}.ParaMap_DCEmask_align.Data(YRange(v,:),XRange(v,:),ZRange(v,:),:));
    
    roi      = double(ROI_multiV{v}.ParaMap_DCEmask_align.Data(YRange(v,:),XRange(v,:),ZRange(v,:)));
    roi(roi>1) = 1;
    roi(isnan(roi)) = 0;
    roi(isinf(roi)) = 0;
    se = strel('disk',2);
    roi = imdilate(roi,se);
    roi = imfill(roi,'holes');

    %adc = ADC_multiV{v}.ParaMap_ADC_reg.Data(YRange(v,:),XRange(v,:),ZRange(v,:));
    adc = ADC_multiV{v}.Calculated_ADC_reg(YRange(v,:),XRange(v,:),ZRange(v,:));
    adc(isnan(adc)) = 0;
    adc(isinf(adc)) = 0;
    
    
    if v == Visit_num_f
        savename = [subname 'T',num2str(Visit_num_f-1) '_FixedImagesForReg.mat'];
    else
        savename = [subname 'T',num2str(v-1),'_MovingToBeRegTo_T',num2str(Visit_num_f-1) '.mat'];
    end
    
    save([outpath savename],'avgdce','anatomical','roi','roi_orig','adc','-v7.3')
    
end


%% Convert hi res and tumor mask files to .mhd files and also added blank slice ends to pad the registration to avoid cutting off ends
for v = Visits
    if v == Visit_num_f
        savename = [subname 'T',num2str(Visit_num_f-1) '_FixedImagesForReg.mat'];
    else
        savename = [subname 'T',num2str(v-1),'_MovingToBeRegTo_T',num2str(Visit_num_f-1) '.mat'];
    end
    load([outpath savename])
    
    %pad ends
    mask  = zeros(size(roi,1),size(roi,2),size(roi,3)+2);
    image = zeros(size(roi,1),size(roi,2),size(roi,3)+2);
    mask(:,:, 2:end-1) = roi;
    
    image(:,:,2:end-1) = anatomical; %avgdce;
    image([1,end],:,:) = 0;
    image(:,[1,end],:) = 0;
    
    init = 5*max(image(:))*mask + image; %range 5-10
    init = 255 * init / max(init(:));
    
    %weight range 0.5-5
    image = image .* (1 + weight*imgaussfilt(double(mask), 2));
    image = 255 * image / max(image(:));
    
    % convert images to .mhd files
    Convert_MAT2MHD_ForElastix(mask,  outpath, 'outputname', [subname 'T',num2str(v-1), '_TumorMask_hs'])
    Convert_MAT2MHD_ForElastix(image, outpath, 'outputname', [subname 'T',num2str(v-1), '_hs'])
    Convert_MAT2MHD_ForElastix(init,  outpath, 'outputname', [subname 'T',num2str(v-1), '_Init_hs'])

% %     if ~exist([outpath,subname,'T0regtoT1/'])
% %         mkdir([outpath,subname,'T0regtoT1/']);
% %     end
% %     if ~exist([outpath,subname,'T2regtoT1/'])
% %         mkdir([outpath,subname,'T2regtoT1/']);
% %     end
% %     outpath0 = [outpath,subname,'T0regtoT1/'];
% %     outpath2 = [outpath,subname,'T2regtoT1/'];

% %     % convert images to .mhd files
% %     Convert_MAT2MHD_ForElastix(mask,  outpath0, 'outputname', [subname 'T',num2str(v-1), '_TumorMask_hs'])
% %     Convert_MAT2MHD_ForElastix(image, outpath0, 'outputname', [subname 'T',num2str(v-1), '_hs'])
% %     Convert_MAT2MHD_ForElastix(init,  outpath0, 'outputname', [subname 'T',num2str(v-1), '_Init_hs'])
% % 
% %     Convert_MAT2MHD_ForElastix(mask,  outpath2, 'outputname', [subname 'T',num2str(v-1), '_TumorMask_hs'])
% %     Convert_MAT2MHD_ForElastix(image, outpath2, 'outputname', [subname 'T',num2str(v-1), '_hs'])
% %     Convert_MAT2MHD_ForElastix(init,  outpath2, 'outputname', [subname 'T',num2str(v-1), '_Init_hs'])

end
