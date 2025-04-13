% clear
% 
% load('/Volumes/ChengyueWu2/DataBackUp/OMG_Assist/ClinicalPipeline/Data/Processed/337_b3_DCE_align.mat')
% load('/Volumes/ChengyueWu2/DataBackUp/OMG_Assist/ClinicalPipeline/Data/Processed/337_b3_DWI_align.mat')
% load('/Volumes/ChengyueWu2/DataBackUp/OMG_Assist/ClinicalPipeline/Data/Processed/337_b3_LocInfo.mat')
% 
% 
% % function [] = RigidRegister_2D_DWI()
% 
% esize = 0; 
% DCE_Data = DCE_align.Data(:,:,end:(-1):1,1);
% Refer_org= DCE_Data;
% % 
% Target_org = DWI_align.Data(:,:,:,1); %135:270,85:200
% Target_org(Target_org < 0)    = 0;
% Target_org(isnan(Target_org)) = 0;
% 
% 
% Target_Region = '50:end,60:190,slc_t1m';
% zidx_valid_1   = zidx_DWI2Dynamic_1;
% zidx_valid_end = zidx_DWI2Dynamic_end;
% xshift = 50 - 1;
% yshift = 60 - 1;
% XWorldLimits = [1, xsize] - 0.5 - xshift;
% YWorldLimits = [1, ysize] - 0.5 - yshift;
% RA_whole = imref2d([ysize, xsize]); %, XWorldLimits, YWorldLimits
% 
% % % % Target_org = T1Mapping_align.Data(:,:,:,1); %135:270,85:200
% % % % Target_org(Target_org < 0)    = 0;
% % % % Target_org(isnan(Target_org)) = 0;
% 
% % zidx_valid_1  = 1;
% % zidx_valid_end= zsize;
% % zidx_valid_1   = max(1,     zidx_DTI2Dynamic_1  -esize);
% % zidx_valid_end = min(zsize, zidx_DTI2Dynamic_end+esize);
% 
% eval(['Target_org_fcs = Target_org(',Target_Region,');'])
% eval(['Refer_org_fcs  = Refer_org (',Target_Region,');'])
% RA_fcs = imref2d(size(Refer_org_fcs));
% 
% [optimizer, metric] = imregconfig('multimodal');
% metric.NumberOfSpatialSamples = 1000;
% optimizer.InitialRadius = 0.0003;
% optimizer.GrowthFactor  = 1.0001;
% optimizer.MaximumIterations = 500;
%       
% transform       = imregtform(Target_org_fcs, Refer_org_fcs, 'rigid', optimizer, metric);
% Target_reg_fcs  = imwarp(Target_org_fcs, transform);
% Target_reg_fcs2 = imwarp(Target_org_fcs, transform,'OutputView',RA_fcs);
% Target_reg_slice= imwarp(Target_org(:,:,slc_t1m), transform,'OutputView',RA_whole);
% 
% 
% Target_reg = zeros(size(Target_org));
% for i = zidx_valid_1:zidx_valid_end
%     Target_reg_slice = imwarp(Target_org(:,:,i), transform,'OutputView',RA_slice);
%     Target_reg(:,:,i)= Target_reg_slice;
% end
% % end
% 
% %% 
% for i = slc_t1m %:zidx_valid_end
%     subplot(1,2,1)
%     imshowpair(Target_reg(:,:,i), Refer_org(:,:,i))
%     subplot(1,2,2)
%     imshowpair(Target_reg_slice, Refer_org(:,:,i)) %Target_reg(:,:,i)
%     title(i)
%     pause
% end

%%

function [Target_reg_slice, transform2D, RA_full] = RigidRegister_2D(Target_org_slice, Refer_org_slice, TargetRegion_YX)
    
% Target_org  = zeros(256,256);
% Refer_org = zeros(256,256);
% Target_org(50:end, 60:190) = DWI_align.Data(50:end, 60:190,slc_t1m,1);
% Refer_org(50:end, 60:190)  = DCE_Data(50:end, 60:190,slc_t1m);

[ysize, xsize] = size(Target_org_slice);
if isempty(TargetRegion_YX)
    TargetRegion_YX = [1, ysize; 1, xsize];
end

Target_org_fcs = Target_org_slice(TargetRegion_YX(1,1):TargetRegion_YX(1,2),...
                                  TargetRegion_YX(2,1):TargetRegion_YX(2,2));
Refer_org_fcs= Refer_org_slice(   TargetRegion_YX(1,1):TargetRegion_YX(1,2),...
                                  TargetRegion_YX(2,1):TargetRegion_YX(2,2));
CenterY = round(sum(TargetRegion_YX(1,:)) / 2);
CenterX = round(sum(TargetRegion_YX(2,:)) / 2);
YWLimit_fcs = TargetRegion_YX(1,:) - CenterY;
XWLimit_fcs = TargetRegion_YX(2,:) - CenterX;

YWLimit_full= [1, ysize]- CenterY;
XWLimit_full= [1, xsize]- CenterX;

RA_fcs = imref2d(size(Refer_org_fcs), XWLimit_fcs, YWLimit_fcs);
RA_full= imref2d([ysize, xsize], XWLimit_full, YWLimit_full);

[optimizer, metric] = imregconfig('multimodal');
metric.NumberOfSpatialSamples = 1000;
optimizer.InitialRadius = 0.0003;
optimizer.GrowthFactor  = 1.0001;
optimizer.MaximumIterations = 500;
      
transform2D     = imregtform(Target_org_fcs, RA_fcs, Refer_org_fcs, RA_fcs, 'rigid', optimizer, metric);
Test_reg_fcs    = imwarp(Target_org_fcs,   RA_fcs,  transform2D,'OutputView',RA_fcs);
Target_reg_slice= imwarp(Target_org_slice, RA_full, transform2D,'OutputView',RA_full);

%%
figure(1)
imshowpair(Target_org_slice, Refer_org_slice)

figure(2)
imshowpair(Test_reg_fcs, Refer_org_fcs)
title('focus on the interested region')

figure(3)
imshowpair(Target_reg_slice(TargetRegion_YX(1,1):TargetRegion_YX(1,2), TargetRegion_YX(2,1):TargetRegion_YX(2,2)),...
           Refer_org_slice( TargetRegion_YX(1,1):TargetRegion_YX(1,2), TargetRegion_YX(2,1):TargetRegion_YX(2,2)))
figure(4)
imshowpair(Target_reg_slice,Refer_org_slice)

drawnow



