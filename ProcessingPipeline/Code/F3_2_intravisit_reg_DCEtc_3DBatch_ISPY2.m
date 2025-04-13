%% example 2: Intra-visit alignment - DCE-time course registration

clear
clc
startpath = pwd;

% % matpath = '/net/bananastand/data1/MDACC_breast/Processed/';
matpath    = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';

Patient_list = {
'104268'...
};
%%%%% for convinience of post-processing, all formated data are moved to a
%%%%% single folder

for patientID = 1:numel(Patient_list)
    PatientName = Patient_list{patientID};
    datapath    = [matpath,   PatientName,  '/'];
        
    figpath  = [matpath, PatientName, '/figures/'];
    if ~isfolder(figpath)
        mkdir(figpath)
    end
    
    for scanID = 0 %0:2 
        CaseNum = [PatientName 'T' num2str(scanID)]
    
        %%
        load([datapath CaseNum,'_DCE_align.mat'])
        load([datapath CaseNum,'_TargetRegion_auto.mat'])

        %first frame of DCE as reference
        DCE_data = DCE_align.Data; 
        Refer_org= DCE_data(:,:,:,1);

        res_YXZ = [DCE_align.resolution; DCE_align.slicethickness];

        %% Registration over DCE time course
        NumT = size(DCE_data,4);
        tform3D_DCEtc_cell = cell(1, NumT);
        RA_DCEtc_cell      = cell(1, NumT);

        DCE_reg   = DCE_align;
        for t = 2:NumT
            t
            Target_org = DCE_data(:,:,:,t); 
            Target_org(Target_org < 0)    = 0;
            Target_org(isnan(Target_org)) = 0;

            [Target_reg, tform3D_dce, RA_dce] = RigidRegister_3D(Target_org, Refer_org, TargetRegion_YXZ, res_YXZ);
            DCE_reg.Data(:,:,:,t) = Target_reg;

            tform3D_DCEtc_cell{t} = tform3D_dce;
            RA_DCEtc_cell{t}      = RA_dce;

        end
        dce_data = DCE_reg.Data;
        dce_data(dce_data > ((2^15)-1)) = 2^15-1;
        dce_data(dce_data < 0) = 0;
        DCE_reg.Data = int16(dce_data);

            %% %%% plot (optional)
            %
            close all
            h1 = figure(1);
            set(h1, 'Unit','characters','Position', [10, 30, 250, 50]);

            z = round(mean(TargetRegion_YXZ(3,:)));
            for t = 2:NumT 
                subplot(1,2,1)
                imshowpair(DCE_align.Data(:,:,z,t), DCE_align.Data(:,:,z,1))
                title('original')
                ax = gca;
                ax.FontSize = 22;
                ax.FontName = 'Times New Roman';

                subplot(1,2,2)
                imshowpair(DCE_reg.Data(:,:,z,t),  DCE_reg.Data(:,:,z,1))
                title('registered')
                ax = gca;
                ax.FontSize = 22;
                ax.FontName = 'Times New Roman';

                saveas(h1, [figpath CaseNum, '_DCEtc_3DRigidReg_t',num2str(t),'.jpg'])
            end
            %}

        %%
        save([datapath, CaseNum,'_DCEtc_3DRigidReg.mat'], 'DCE_reg','RA_DCEtc_cell','tform3D_DCEtc_cell','-v7.3')

    end
end

beep


