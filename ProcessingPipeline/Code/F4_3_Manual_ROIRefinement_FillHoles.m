clc; clearvars; close all
startpath = pwd;

matpath    = '/Volumes/SanDiskSSD/ISPY2/Data/Processed/';
allpath = '/Users/rp33847/UserFolders/ISPY2/ByPatient/';

%{
%}

Patient_list = {
'698094'...
};

for pidx = 1:length(Patient_list)

    PatientName = Patient_list{pidx};
    disp(PatientName); disp('   ')

    savepath = [allpath PatientName '/'];
    unregpath  = [matpath,   PatientName,  '/UnregisteredData/'];
    final = 0;
    
    for scanID = 0
        while final==0
            close all
            CaseNum = [PatientName 'T' num2str(scanID)];
    
            load([unregpath CaseNum,'_DCEtc_3DRigidReg.mat'],'DCE_reg')
            load([unregpath CaseNum '_DCEmask_Manual_FCM_3DRigidReg.mat'], 'TumorROI','TumorMargin','U','V','threshold','LMM')
            load([unregpath CaseNum,'_ParaMap_DWImask_3DRigidReg.mat'])
    
            load([unregpath CaseNum,'_TargetRegion_auto.mat'])
            load([unregpath CaseNum '_LocInfo.mat'])
    
            %% Fix Tumor ROI
            slices1 = find(sum(sum(TumorROI)))';
            slices = [slices1(1)-1 slices1 slices1(end)+1];
            figure(1)
            imagesc(imtile(TumorROI(:,:,slices))); axis image
    
            figure(2)
            imagesc(imtile(DCE_reg.Data(:,:,slices,2))); axis image
    
            figure(3)
            combined = double(DCE_reg.Data(:,:,:,2)) + double(TumorROI).*(max(double(DCE_reg.Data),[],'all'));
            imagesc(imtile(combined(:,:,slices))); colorbar; axis image
    
            figure(4)
            ROI = zeros(size(TumorROI));
            for slc=slices
                imagesc(imtile(combined(:,:,slc))); colorbar; clim([0 0.9*max(DCE_reg.Data,[],'all')]); axis image
                title([num2str(slc) ' / ' num2str(slices(end))]);
                [ROI(:,:,slc)] = roipoly;
            end
            
            TumorROI_New = TumorROI;
            TumorROI_New = TumorROI_New + ROI;
            for z = 1:size(TumorROI_New,3)
                TumorROI_New(:,:,z) = imfill(TumorROI_New(:,:,z),'holes');
            end
            TumorROI_New(TumorROI_New>=1)=1;
            
            %% Display new tumor ROI
    
            figure(40)
            combined_new = double(DCE_reg.Data(:,:,:,2)) + double(TumorROI_New).*(max(double(DCE_reg.Data),[],'all'));
            imagesc(imtile(combined_new(:,:,slices))); colorbar; axis image

            figure(41)
            testvol_DCE = TumorROI_New;
            testvol_DWI = ParaMap_DWImask_reg.Data;
            
            col_dce=[0.8 0 0];
            col_dwi=[0 0.8 0];
    
            hiso1 = patch(isosurface(testvol_DCE,0),'FaceColor',col_dce,'EdgeColor','none');
            hold on
            hiso2 = patch(isosurface(testvol_DWI,0),'FaceColor',col_dwi,'EdgeColor','none');
            legend({'DCE-derived tumor', 'DWI-derived tumor'})
    
            axis on; grid on
            lighting phong;
            isonormals(testvol_DCE,hiso1);
            isonormals(testvol_DWI,hiso2);
            
            alpha(0.5);
            set(gca,'DataAspectRatio',[1/DCE_reg.resolution(1) 1/DCE_reg.resolution(2) 1/DCE_reg.slicethickness])
            camlight;
            view([0.5 0.8 0.5])
            
            ax = gca;
            ax.FontSize = 16;
            
            set(gcf, 'Unit','characters','Position', [10, 30, 50, 30]);
    
            disp('  ')
            
            final = input('1 for save\n');
        end
        TumorROI = TumorROI_New;
        save([unregpath CaseNum '_DCEmask_Manual_FCM_3DRigidReg.mat'], 'TumorROI','TumorMargin','U','V','threshold','LMM')
        final = 1;
    end  
end
