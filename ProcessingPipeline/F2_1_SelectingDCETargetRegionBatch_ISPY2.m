clear
clc
startpath = pwd;

% % matpath = '/net/bananastand/data1/MDACC_breast/Processed/';
matpath    = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';

VariableNames = {'DWI','DCE','ParaMap_ADC','ParaMap_DWImask'};

Patient_list = {
'104268'...
};

%% %%% for convinience of post-processing, all formated data are moved to a
%%%%% single folder

for patientID = 1:numel(Patient_list)
    PatientName = Patient_list{patientID};
    datapath  = [matpath,   PatientName,  '/'];
        
    for scanID = 0:2 %only have segmentations for first 3 visits
        
        CaseNum = [PatientName 'T' num2str(scanID)]
        
        load([datapath CaseNum,'_ParaMap_DWImask_align.mat'])
        load([datapath CaseNum,'_ParaMap_DCEmask_align.mat'])
        % load([datapath CaseNum '_DCEmask_Manual_FCM_3DRigidReg_Old.mat'], 'TumorROI')

        %Define ROIs and define bounding box
        ROI1 = double(ParaMap_DWImask_align.Data);
        ROI2 = double(ParaMap_DCEmask_align.Data);

        ROI = ROI1 + ROI2;
        %ROI = ROI1;

        ROI(ROI>1) = 1;
        ROI = imfill(ROI,'holes');
        box = regionprops3(ROI,'BoundingBox');
        TargetRegion_YXZ = [floor(box.BoundingBox(2)), floor(box.BoundingBox(2))+ceil(box.BoundingBox(5)); ...
                            floor(box.BoundingBox(1)), floor(box.BoundingBox(1))+ceil(box.BoundingBox(4)); ...
                            floor(box.BoundingBox(3)), floor(box.BoundingBox(3))+ceil(box.BoundingBox(end))];
        
        %
        R = 10; Rz = 5; %%% 10 5 Could make this smaller
        TargetRegion_YXZ(1,1) = max([1, TargetRegion_YXZ(1,1)-R]);  TargetRegion_YXZ(1,2) = min([size(ROI,1), TargetRegion_YXZ(1,2)+R]); 
        TargetRegion_YXZ(2,1) = max([1, TargetRegion_YXZ(2,1)-R]);  TargetRegion_YXZ(2,2) = min([size(ROI,2), TargetRegion_YXZ(2,2)+R]); 
        TargetRegion_YXZ(3,1) = max([1, TargetRegion_YXZ(3,1)-Rz]); TargetRegion_YXZ(3,2) = min([size(ROI,3), TargetRegion_YXZ(3,2)+Rz]); 
        %{
        translate = 0;
        if translate==1
            TargetRegion_YXZ(1,1) = max([1, TargetRegion_YXZ(1,1)+R]);  TargetRegion_YXZ(1,2) = min([size(ROI,1), TargetRegion_YXZ(1,2)+R]); 
            TargetRegion_YXZ(2,1) = max([1, TargetRegion_YXZ(2,1)]);  TargetRegion_YXZ(2,2) = min([size(ROI,2), TargetRegion_YXZ(2,2)]); 
            TargetRegion_YXZ(3,1) = max([1, TargetRegion_YXZ(3,1)-Rz]); TargetRegion_YXZ(3,2) = min([size(ROI,3), TargetRegion_YXZ(3,2)+Rz]);
        elseif translate==2
            TargetRegion_YXZ(1,1) = max([1, TargetRegion_YXZ(1,1)-R]);  TargetRegion_YXZ(1,2) = min([size(ROI,1), TargetRegion_YXZ(1,2)-R]); 
            TargetRegion_YXZ(2,1) = max([1, TargetRegion_YXZ(2,1)]);  TargetRegion_YXZ(2,2) = min([size(ROI,2), TargetRegion_YXZ(2,2)]); 
            TargetRegion_YXZ(3,1) = max([1, TargetRegion_YXZ(3,1)-Rz]); TargetRegion_YXZ(3,2) = min([size(ROI,3), TargetRegion_YXZ(3,2)+Rz]);
        elseif translate==3
            TargetRegion_YXZ(1,1) = max([1, TargetRegion_YXZ(1,1)]);  TargetRegion_YXZ(1,2) = min([size(ROI,1), TargetRegion_YXZ(1,2)]); 
            TargetRegion_YXZ(2,1) = max([1, TargetRegion_YXZ(2,1)+R]);  TargetRegion_YXZ(2,2) = min([size(ROI,2), TargetRegion_YXZ(2,2)+R]); 
            TargetRegion_YXZ(3,1) = max([1, TargetRegion_YXZ(3,1)-Rz]); TargetRegion_YXZ(3,2) = min([size(ROI,3), TargetRegion_YXZ(3,2)+Rz]);
        elseif translate==4
            TargetRegion_YXZ(1,1) = max([1, TargetRegion_YXZ(1,1)]);  TargetRegion_YXZ(1,2) = min([size(ROI,1), TargetRegion_YXZ(1,2)]); 
            TargetRegion_YXZ(2,1) = max([1, TargetRegion_YXZ(2,1)-R]);  TargetRegion_YXZ(2,2) = min([size(ROI,2), TargetRegion_YXZ(2,2)-R]); 
            TargetRegion_YXZ(3,1) = max([1, TargetRegion_YXZ(3,1)-Rz]); TargetRegion_YXZ(3,2) = min([size(ROI,3), TargetRegion_YXZ(3,2)+Rz]);
        elseif translate==5
            TargetRegion_YXZ(1,1) = max([1, TargetRegion_YXZ(1,1)-R]);  TargetRegion_YXZ(1,2) = min([size(ROI,1), TargetRegion_YXZ(1,2)+R]); 
            TargetRegion_YXZ(2,1) = max([1, TargetRegion_YXZ(2,1)-R]);  TargetRegion_YXZ(2,2) = min([size(ROI,2), TargetRegion_YXZ(2,2)+R]); 
            TargetRegion_YXZ(3,1) = max([1, TargetRegion_YXZ(3,1)+Rz]); TargetRegion_YXZ(3,2) = min([size(ROI,3), TargetRegion_YXZ(3,2)+Rz]);
        elseif translate==6
            TargetRegion_YXZ(1,1) = max([1, TargetRegion_YXZ(1,1)-R]);  TargetRegion_YXZ(1,2) = min([size(ROI,1), TargetRegion_YXZ(1,2)+R]); 
            TargetRegion_YXZ(2,1) = max([1, TargetRegion_YXZ(2,1)-R]);  TargetRegion_YXZ(2,2) = min([size(ROI,2), TargetRegion_YXZ(2,2)+R]); 
            TargetRegion_YXZ(3,1) = max([1, TargetRegion_YXZ(3,1)-Rz]); TargetRegion_YXZ(3,2) = min([size(ROI,3), TargetRegion_YXZ(3,2)-Rz]);
        else
            TargetRegion_YXZ(1,1) = max([1, TargetRegion_YXZ(1,1)-R]);  TargetRegion_YXZ(1,2) = min([size(ROI,1), TargetRegion_YXZ(1,2)+R]); 
            TargetRegion_YXZ(2,1) = max([1, TargetRegion_YXZ(2,1)-R]);  TargetRegion_YXZ(2,2) = min([size(ROI,2), TargetRegion_YXZ(2,2)+R]); 
            TargetRegion_YXZ(3,1) = max([1, TargetRegion_YXZ(3,1)-Rz]); TargetRegion_YXZ(3,2) = min([size(ROI,3), TargetRegion_YXZ(3,2)+Rz]);
        end
        
        %}

        TargetRegion_YXZ
        save([datapath CaseNum '_TargetRegion_auto.mat'], 'TargetRegion_YXZ')

    end
end
%}
