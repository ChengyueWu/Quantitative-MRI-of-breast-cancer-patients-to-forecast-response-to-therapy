clc; clearvars; close all
startpath = pwd;

matpath    = '/Volumes/SanDiskSSD/ISPY2/Data/Processed/';
allpath = '/Users/rp33847/UserFolders/ISPY2/ByPatient/';

%{
'928170'...
'963303'...
'973837'...

%}

Patient_list = {
'698094'...
};
thresh = 0.5;

for pidx = 1:length(Patient_list)
disp('Outline sections to remove')
    PatientName = Patient_list{pidx};
    disp(PatientName); disp('   ')

    savepath = [allpath PatientName '/'];

    unregpath  = [matpath,   PatientName,  '/UnregisteredData/'];
    
    final = 0;
    for scanID = 2
        while final==0
            close all
            CaseNum = [PatientName 'T' num2str(scanID)];
    
            load([unregpath CaseNum,'_DCEtc_3DRigidReg.mat'],'DCE_reg')
            load([unregpath CaseNum '_DCEmask_Manual_FCM_3DRigidReg.mat'], 'TumorROI','TumorMargin','U','V','threshold','LMM')
            load([unregpath CaseNum,'_ParaMap_DWImask_3DRigidReg.mat'])
    
            load([unregpath CaseNum,'_TargetRegion_auto.mat'])
            load([unregpath CaseNum '_LocInfo.mat'])
    
            %% Fix Tumor ROI
            slices = find(sum(sum(TumorROI)))';
            figure(1)
            imagesc(imtile(TumorROI(:,:,slices))); axis image
    
            % figure(2)
            % imagesc(imtile(DCE_reg.Data(:,:,slices,4))); axis image
    
            figure(2)
            combined = double(DCE_reg.Data(:,:,:,4)) + double(TumorROI).*thresh*1e3;
            imagesc(imtile(combined(:,:,slices))); colorbar; clim([0 0.9*max(DCE_reg.Data,[],'all')]); axis image
    
            figure(3)
            ROI = zeros(size(TumorROI));
            for slc=slices
                imagesc(imtile(combined(:,:,slc))); colorbar; clim([0 0.9*max(DCE_reg.Data,[],'all')]); axis image
                title([num2str(slc) ' / ' num2str(slices(end))]);
                [ROI(:,:,slc),~,~] = roipoly;
            end
            
            TumorROI_New = TumorROI;
            TumorROI_New = TumorROI_New.*(1-ROI); %%% CHANGE HERE
            
            %% Display new tumor ROI
    
            figure(4)
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

            figure(2)
            slices = find(sum(sum(TumorROI_New)))';
            Mask   = double(TumorROI_New);
            img    = double(DCE_reg.Data(:,:,:,2));
            %img    = double(DCE_reg.Data(:,:,:,2)-DCE_reg.Data(:,:,:,1));
    
            BW = zeros(size(Mask));
            for z = slices
                bw_z = edge(Mask(:,:,z));
                BW(:,:,z) = bw_z;
            end
            img_disp = repmat(permute(img(:,:,slices) / prctile(img(:),99.9), [1 2 4 3]), [1 1 3 1]);
            img_disp = img_disp .* ~repmat(permute(BW(:,:,slices),[1 2 4 3]), [1,1,3,1]);
            img_disp(:,:,1,:) = img_disp(:,:,1,:) + permute(BW(:,:,slices),[1 2 4 3]);
            montage(img_disp)
            set(gcf,'Units','characters','Position',[10 10 150 80]);
            title(['Tumor ROI ',PatientName,' Visit ',num2str(scanID+1)],'FontSize',25);
    
            disp('  ')
            
            final = input('1 for save\n');
        end
        TumorROI = TumorROI_New;
        save([unregpath CaseNum '_DCEmask_Manual_FCM_3DRigidReg.mat'], 'TumorROI','TumorMargin','U','V','threshold','LMM')
    end  
end
