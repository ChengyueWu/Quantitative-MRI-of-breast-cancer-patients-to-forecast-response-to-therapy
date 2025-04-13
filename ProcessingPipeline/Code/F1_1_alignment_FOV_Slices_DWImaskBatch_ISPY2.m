clear
clc
startpath = pwd;

% % matpath = '/net/bananastand/data1/MDACC_breast/Processed/';
matpath    = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';

Patient_list = {
'104268'...
};

%% %%% for convinience of post-processing, all formated data are moved to a single folder

for patientID = 1:numel(Patient_list)
    PatientName = Patient_list{patientID};
    
    datapath = [matpath, PatientName, '/']; % matpath; %
    ROIpath  = datapath; 
    
    figpath  = [datapath, 'figures/'];
    if ~isfolder(figpath)
        mkdir(figpath)
    end
    
    for scanID = 0:2 %only have segmentations for first 3 visits
        
        CaseNum = [PatientName 'T' num2str(scanID)]

        %% load DWI and DWI-derived mask
        load([datapath CaseNum, '_DWI_scandata.mat'])
        load([datapath CaseNum, '_ParaMap_DWImask_scandata.mat'])

        Zlist_DWI_whole = DWI.pos(:,3,1);
        Zlist_DWI_mask  = ParaMap_DWImask.pos(:,3);

        %% find slice indices of mask
        ZIDX_mask2whole = [];
        for z_idx = 1:length(Zlist_DWI_mask)
            [~,m2w] = min(abs(Zlist_DWI_whole - Zlist_DWI_mask(z_idx)));
            ZIDX_mask2whole = [ZIDX_mask2whole; m2w];
        end
        
        %% convert mask to the DWI grid
        ParaMap_DWImask_whole = DWI;
        ParaMap_DWImask_whole.Data = zeros(size(DWI.Data(:,:,:,1)));
        for z_idx = 1:length(Zlist_DWI_mask)
            ParaMap_DWImask_whole.Data(:,:,ZIDX_mask2whole(z_idx)) = ParaMap_DWImask.Data(:,:,z_idx);
        end
        ParaMap_DWImask_whole.pos  = ParaMap_DWImask_whole.pos(:,:,1);
        ParaMap_DWImask_whole.ori  = ParaMap_DWImask_whole.ori(:,:,1);
        ParaMap_DWImask_whole.slloc= ParaMap_DWImask_whole.slloc(:,:,1);
        ParaMap_DWImask_whole.dicom_header = ParaMap_DWImask.dicom_header;
        
            %% visualize
            close all
            slices = ZIDX_mask2whole';
            Mask   = ParaMap_DWImask_whole.Data;
            img    = DWI.Data(:,:,:,2);
    
            BW = zeros(size(Mask));
            for z = slices
                bw_z = edge(Mask(:,:,z));
                BW(:,:,z) = bw_z;
            end
            img_disp = repmat(permute(img(:,:,slices) / prctile(img(:),99.9), [1 2 4 3]), [1 1 3 1]);
            img_disp = img_disp .* ~repmat(permute(BW(:,:,slices),[1 2 4 3]), [1,1,3,1]);
            img_disp(:,:,1,:) = img_disp(:,:,1,:) + permute(BW(:,:,slices),[1 2 4 3]);
            figure(1);
            montage(img_disp)
            set(gcf,'Units','characters','Position',[10 10 150 80]);
            
            pause(1)
            figname = [CaseNum, '_Montage_DWImask_whole'];
            saveas(gcf, [figpath, figname, '.jpg'])

            
        %% save
        savename = [CaseNum, '_ParaMap_DWImask_whole.mat'];
        save([datapath savename], 'ParaMap_DWImask_whole')
    end
end
