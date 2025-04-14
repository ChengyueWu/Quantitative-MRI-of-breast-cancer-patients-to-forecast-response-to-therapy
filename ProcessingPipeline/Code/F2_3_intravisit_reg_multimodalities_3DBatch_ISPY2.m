%% example 2: Intra-visit alignment - multi-modal registration (focus on tumor area)

clear
clc
startpath = pwd;

matpath    = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';

fun_h = @(para, xdata)(xdata>=para(3)) .* (2*para(1)./(1 + exp(-para(2).*(xdata-para(3)))) - para(1));
option = optimoptions('lsqcurvefit', 'Display', 'off');

Patient_list = {
'104268'...
};


%% %%% for convinience of post-processing, all formated data are moved to a
for patientID = 1:numel(Patient_list)
    PatientName = Patient_list{patientID};
    datapath  = [matpath, PatientName,  '/'];
        
    figpath  = [matpath, PatientName, '/figures/'];
    if ~isfolder(figpath)
        mkdir(figpath)
    end
    
    for scanID = 1 %0:2
        CaseNum = [PatientName 'T' num2str(scanID)]
        
        clear DCE_align DWI_align ParaMap_ADC_align ParaMap_DCEmask_align ParaMap_DWImask_align
        clear DCE_reg   DWI_reg   ParaMap_ADC_reg   ParaMap_DCEmask_reg   ParaMap_DWImask_reg
        
        %%
        load([datapath CaseNum,'_DCEtc_3DRigidReg.mat'],'DCE_reg')
        load([datapath CaseNum,'_ParaMap_DCEmask_align.mat'])
        
        load([datapath CaseNum,'_DWI_align.mat'])
        load([datapath CaseNum,'_ParaMap_DWImask_align.mat'])

        load([datapath CaseNum,'_TargetRegion_auto.mat'])

        
        %% first frame of DCE as reference
        DCE_data = double(DCE_reg.Data); 
        
        TotalInt = squeeze(sum(sum(sum(DCE_data .* repmat(ParaMap_DWImask_align.Data, [1,1,1,size(DCE_data,4)])))));
        TotalInt = TotalInt - TotalInt(1);
        para0 = [max(TotalInt), 0.5, 5];

        %%%% fitting half-logistic
        para_h = lsqcurvefit(fun_h, para0, 1:length(TotalInt), TotalInt', [0,0,1], [2*max(TotalInt),2*max(TotalInt),length(TotalInt)], option);
        estT0  = floor(max([1, para_h(3)]));  
        EnhT = min([estT0+1, size(DCE_reg.Data)]);
            
        res_YXZ = [DCE_reg.resolution; DCE_reg.slicethickness];

        
        %% Registration for DWI
        %
        DWI_Data = double(DWI_align.Data);
        DWI_Data(isnan(DWI_Data)) = 0;
        intensities = [sum(sum(sum(DWI_Data(:,:,:,1)))) sum(sum(sum(DWI_Data(:,:,:,2)))) ...
            sum(sum(sum(DWI_Data(:,:,:,3)))) sum(sum(sum(DWI_Data(:,:,:,4))))];
        [intensities_sorted,I] = sort(intensities,'descend');
        bval100 = find(I==2);
        dwi_org = DWI_align.Data(:,:,:,bval100); 
        dwi_org(dwi_org < 0)    = 0;
        dwi_org(isnan(dwi_org)) = 0;

        Refer_org= DCE_data(:,:,:,2);
        
        %% New!
        dwi_org  = double(ParaMap_DWImask_align.Data);
        Refer_org= double(ParaMap_DCEmask_align.Data);

        [~, tform3D_DWI, RA_dwi] = RigidRegister_3D(dwi_org, Refer_org, TargetRegion_YXZ, res_YXZ);

        DWI_reg = DWI_align;
        for f = 1:size(DWI_reg.Data,4)
            Target_reg = imwarp(DWI_align.Data(:,:,:,f), RA_dwi, tform3D_DWI,'OutputView',RA_dwi);
            DWI_reg.Data(:,:,:,f)= Target_reg;
        end
        dwi_data = DWI_reg.Data;
        dwi_data(dwi_data > ((2^15)-1)) = 2^15-1;
        dwi_data(dwi_data < 0) = 0;
        DWI_reg.Data = int16(dwi_data);

            %%%%% plot (optional)
            %
            h1 = figure(1);
            set(h1, 'Unit','characters','Position', [10, 30, 250, 50]);
            
            z = 55;
            for f = 1:size(DWI_align.Data,4)
                figure(1)
                subplot(1,2,1)
                imshowpair(double(DWI_align.Data(:,:,z,f)), ...
                           double(DCE_reg.Data(:,:,z,EnhT)))
                ax = gca;
                ax.FontSize = 22;

                title('original')
                subplot(1,2,2)
                imshowpair(double(DWI_reg.Data(:,:,z,f)), ...
                           double(DCE_reg.Data(:,:,z,EnhT)))
                title('registered')
                ax = gca;
                ax.FontSize = 22;
                saveas(h1, [figpath CaseNum, '_DWI_3DRigidReg_f',num2str(f),'.jpg'])
            end
            % close all
            %}
        

        Target_reg = imwarp(ParaMap_DWImask_align.Data, RA_dwi, tform3D_DWI,'OutputView',RA_dwi);
        ParaMap_DWImask_reg.Data = Target_reg;
        ParaMap_DWImask_reg.Data = int16(ParaMap_DWImask_reg.Data > 0.5);
                
                %
                h3 = figure(3);
                set(h3, 'Unit','characters','Position', [10, 30, 260, 50]);
                
                z = 55; 
                subplot(1,3,1)
                imshowpair(double(ParaMap_DWImask_align.Data(:,:,z)), ...
                           double(DCE_reg.Data(:,:,z,EnhT)))
                title('original')
                ax = gca;
                ax.FontSize = 22;

                subplot(1,3,2)
                imshowpair(double(ParaMap_DWImask_reg.Data(:,:,z)), ...
                           double(DCE_reg.Data(:,:,z,EnhT)))
                title('registered')
                ax = gca;
                ax.FontSize = 22;

                subplot(1,3,3)
                imshowpair(double(ParaMap_DWImask_reg.Data(:,:,z)), ...
                           double(DCE_reg.Data(:,:,z,1)))
                title('registered (vs t = 0)')
                ax = gca;
                ax.FontSize = 22;

                saveas(h3, [figpath CaseNum, '_ParaMap_DWImask_3DRigidReg.jpg'])
                
                % close all
                %} 

        save([datapath CaseNum '_DWI_3DRigidReg.mat'],'DWI_reg','RA_dwi','tform3D_DWI')
        save([datapath CaseNum '_ParaMap_DWImask_3DRigidReg.mat'],'ParaMap_DWImask_reg','RA_dwi','tform3D_DWI')
        %}

        
    end
end
%end of file

