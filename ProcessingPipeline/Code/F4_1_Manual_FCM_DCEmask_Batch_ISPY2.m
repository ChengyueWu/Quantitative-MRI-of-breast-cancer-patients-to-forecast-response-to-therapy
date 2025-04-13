clear
clc
startpath = pwd;

matpath    = '/Volumes/SanDiskSSD/ISPY2/Data/Processed/';

%VariableNames = {'DWI','DCE','ParaMap_ADC','ParaMap_DWImask'};

%% Patient list
%{
'520471'...
'800447'...
%}

Patient_list = {
'698094'...
};

%%
nonspecific = 0;
 
%slice_list = [];

close all

%%%%% for convinience of post-processing, all formated data are moved to a
%%%%% single folder

for pidx = 1:numel(Patient_list)
    PatientName = Patient_list{pidx};
    outpath  = [matpath,   PatientName,  '/UnregisteredData/'];
    
    figpath  = [matpath, PatientName, '/figures/'];
    if ~isfolder(figpath)
        mkdir(figpath)
    end
    
    for scanID = 2

        CaseNum = [PatientName 'T' num2str(scanID)]
        
        load([outpath CaseNum,'_ParaMap_DWImask_3DRigidReg.mat'])
        TumorROI = double(ParaMap_DWImask_reg.Data);
        load([outpath CaseNum '_DCEmask_Manual_FCM_3DRigidReg.mat'], 'TumorROI')
        
        slices = find(sum(sum(TumorROI)))'; pad = 0;
        slices = 50:55;
        slc_start = min(slices); slc_end = max(slices);
        slices = (slc_start-pad):(slc_end+pad);
        slc1 = (slc_start-pad); slc2 = (slc_end+pad);
        
        clearvars TumorROI
        
        load([outpath CaseNum,'_ParaMap_DWImask_3DRigidReg.mat'])
        load([outpath CaseNum,'_DCEtc_3DRigidReg.mat'])
        
        %% Define the approximate margin of tumor
        %Select slices with tumors
        postcontrast = double(mean(DCE_reg.Data(:,:,:,2:end),4));
        %postcontrast = double(DCE_reg.Data(:,:,:,2));
        precontrast = double(DCE_reg.Data(:,:,:,1));
        img = postcontrast - precontrast;
        figure(100)
        imagesc(imtile(img)); axis image; clim([0 1500])

        figure(101)
        imagesc(imtile(precontrast)); axis image; %clim([0 1500])

        %Generate enhanced image and draw bounding box
        if nonspecific
            postcontrast = double(mean(DCE_reg.Data(:,:,:,2:end),4));
            %postcontrast = double(DCE_reg.Data(:,:,:,2));
            precontrast = double(DCE_reg.Data(:,:,:,1));
            img = postcontrast - precontrast;
            %img_tumor = img(:,:,slc1:slc2);
            img_tumor = img;
            img_tumor = img / prctile(img(:),99);
            img_enh = sum(img_tumor,3);
            
            TumorMargin = zeros(size(img));

            figure(106)
            ex = ceil(mean([slc1 slc2]));
            imagesc(img_tumor(:,:,ex)); axis image
            imagesc(postcontrast(:,:,ex)); axis image
            title([num2str(ex-slc1+1) ' / ' num2str(slc2-slc1+1)],'Fontsize',20)
            [roi,~,~] = roipoly;
            for ii=slc1:slc2
                TumorMargin(:,:,ii) = roi;
            end
        
        else
            postcontrast = double(mean(DCE_reg.Data(:,:,:,2:end),4));
            postcontrast = double(DCE_reg.Data(:,:,:,2));
            precontrast = double(DCE_reg.Data(:,:,:,1));
            img = postcontrast - precontrast;
            %img = postcontrast;
            %img_tumor = img(:,:,slc1:slc2);
            img_tumor = img;
            %img_tumor = img / prctile(img(:),99);
            img_enh = sum(img_tumor,3);
            
            TumorMargin = zeros(size(img));

            figure(106)
            for ii=slc1:slc2
                imagesc(img_tumor(:,:,ii)); axis image
                imagesc(postcontrast(:,:,ii)); axis image
                title([num2str(ii-slc1+1) ' / ' num2str(slc2-slc1+1)],'Fontsize',20)
                %clim([0 2])
                [roi,~,~] = roipoly;
                TumorMargin(:,:,ii) = roi;
            end
    
%             for n = (slc1-3):(slc2+3)
%                 TumorMargin(:,:,n) = roi;
%             end
        end

        %% Enhancement by dividing the post-contrast intensity by the pre-contrast
        [ysize, xsize, zsize, numdyn] = size(DCE_reg.Data);
        %numdyn = 2;

        Ienh = permute(double(DCE_reg.Data(:,:,:,1:numdyn)) .* repmat(TumorMargin, [1,1,1,numdyn]),[1 2 4 3]);
        Ienh = Ienh - repmat(Ienh(:,:,1,:),[1 1 numdyn 1]);
        Ienh = reshape(permute(Ienh,[1 2 4 3]),[ysize*xsize*zsize, numdyn]);

        Condition = sum(Ienh,2);
        Condition(isnan(Condition)) = 0;

        vox = find(Condition);
        X = Ienh(vox,:);

        %% 3. 2-category FCM
        [V,U,~] = FCMClust(X,2);

        %% 4. binarization of lesion membership
        threshold = 0.2; %%% 0.5 CHANGE to 0.2, 0.3
        LMM = zeros(ysize*xsize*zsize,1);
        if sum(V(1,:).^2) > sum(V(2,:).^2)
            LMM(vox) = U(1,:);
        else 
            LMM(vox) = U(2,:);
        end

        LMM = reshape(LMM,[ysize xsize zsize]);

        %% 5. reduce the suprious structure
        BW_lesion = imgaussfilt3(LMM,2)>threshold;
        BW_lesion_ei = bwareaopen(BW_lesion, 60);

        %% 6. hole-filling operation
        BW_lesion_hf = zeros(ysize,xsize,zsize);
        for z = 1:zsize
            BW_lesion_hf(:,:,z) = imfill(BW_lesion_ei(:,:,z),'holes');
        end

        %% final
        TumorROI = BW_lesion_hf;
        %TumorROI = ParaMap_DWImask_reg.Data;
         
        %% visualization
        close all
        
        testvol_DCE = TumorROI;
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
        figname = [CaseNum, '_3DROI_DWI-DCE_Manual'];
        %saveas(gcf, [figpath, figname, '.jpg'])
        
        %%%%%%%%%%%%%%
        figure(2)
        slices = find(sum(sum(TumorROI)))';
        Mask   = double(TumorROI);
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

        figname = [CaseNum, '_Montage_DCEmask_FCM_Manual'];
        %saveas(gcf, [figpath, figname, '.jpg'])

        savename = [PatientName,'_04_V' num2str(scanID+1) '_TumorROIMontage'];
        %saveas(gcf, ['/Users/rp33847/UserFolders/ISPY2/ByPatient/Montage/', savename, '.jpg'])

    %% save
    save([outpath CaseNum '_DCEmask_Manual_FCM_3DRigidReg_Old.mat'], 'TumorROI','TumorMargin','U','V','threshold','LMM')

    end
end
