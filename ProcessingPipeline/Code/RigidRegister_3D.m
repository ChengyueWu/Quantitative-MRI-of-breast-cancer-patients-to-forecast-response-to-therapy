
function [Target_reg, transform3D, RA_full] = RigidRegister_3D(Target_org, Refer_org, TargetRegion_YXZ, res_YXZ)
    
% Target_org  = zeros(256,256);
% Refer_org = zeros(256,256);
% Target_org(50:end, 60:190) = DWI_align.Data(50:end, 60:190,slc_t1m,1);
% Refer_org(50:end, 60:190)  = DCE_Data(50:end, 60:190,slc_t1m);

[ysize, xsize, zsize] = size(Target_org);
if isempty(TargetRegion_YXZ)
    TargetRegion_YXZ = [1, ysize; 1, xsize; 1, zsize];
end

Target_org_fcs = Target_org(TargetRegion_YXZ(1,1):TargetRegion_YXZ(1,2),...
                            TargetRegion_YXZ(2,1):TargetRegion_YXZ(2,2),...
                            TargetRegion_YXZ(3,1):TargetRegion_YXZ(3,2));
                        
Refer_org_fcs  = Refer_org( TargetRegion_YXZ(1,1):TargetRegion_YXZ(1,2),...
                            TargetRegion_YXZ(2,1):TargetRegion_YXZ(2,2),...
                            TargetRegion_YXZ(3,1):TargetRegion_YXZ(3,2));
                        
CenterY = round(sum(TargetRegion_YXZ(1,:)) / 2);
CenterX = round(sum(TargetRegion_YXZ(2,:)) / 2);
CenterZ = round(sum(TargetRegion_YXZ(3,:)) / 2);
YWLimit_fcs = (TargetRegion_YXZ(1,:) - CenterY) * res_YXZ(1);
XWLimit_fcs = (TargetRegion_YXZ(2,:) - CenterX) * res_YXZ(2);
ZWLimit_fcs = (TargetRegion_YXZ(3,:) - CenterZ) * res_YXZ(3);

YWLimit_full= ([1, ysize]- CenterY) * res_YXZ(1);
XWLimit_full= ([1, xsize]- CenterX) * res_YXZ(2);
ZWLimit_full= ([1, zsize]- CenterZ) * res_YXZ(3);

RA_fcs = imref3d(size(Refer_org_fcs),   XWLimit_fcs,  YWLimit_fcs,  ZWLimit_fcs);
RA_full= imref3d([ysize, xsize, zsize], XWLimit_full, YWLimit_full, ZWLimit_full);

[optimizer, metric] = imregconfig('multimodal');
metric.NumberOfSpatialSamples = 5000;
optimizer.InitialRadius = 0.0003;
optimizer.GrowthFactor  = 1.0001;
optimizer.MaximumIterations = 500;
      
transform3D    = imregtform(Target_org_fcs, RA_fcs, Refer_org_fcs, RA_fcs, 'rigid', optimizer, metric);

Target_reg_fcs = imwarp(Target_org_fcs, RA_fcs,  transform3D, 'OutputView',RA_fcs);
Target_reg     = imwarp(Target_org,     RA_full, transform3D, 'OutputView',RA_full);






