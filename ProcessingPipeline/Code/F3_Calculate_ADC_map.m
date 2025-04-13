% calculate parametric maps (for example, ADC) - recommanded to performed after intra-visit
% reg and before inter-visit reg
clearvars; clc
close all

matpath    = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';

Patient_list = {
'104268'...
};

for pidx = 1:length(Patient_list)

    PatientName = Patient_list{pidx};
    datapath  = [matpath,   PatientName,  '/'];  
    figpath   = [matpath, PatientName, '/figures/'];
    
    for scanID = 0:2
        CaseNum = [PatientName 'T' num2str(scanID)]

        load([datapath CaseNum '_DWI_3DRigidReg.mat'],'DWI_reg','RA_dwi','tform3D_DWI')
        %load([outpath CaseNum '_ParaMap_ADC_3DRigidReg.mat'],'ParaMap_ADC_reg')
        %load([outpath CaseNum '_Calculated_ADC_3DRigidReg.mat'],'Calculated_ADC_reg')
        %Provided_ADC_reg = ParaMap_ADC_reg.Data;

        %load([outpath CaseNum '_ParaMap_DWImask_3DRigidReg.mat'],'ParaMap_DWImask_reg')
        %load([outpath CaseNum,'_DCEtc_3DRigidReg.mat'],'DCE_reg')
        %load([outpath CaseNum,'_DWI_trim.mat'])
        %load([outpath CaseNum,'_ParaMap_DWImask_trim.mat'])
        
        load([datapath CaseNum,'_TargetRegion_auto.mat'])
        z = round(mean(TargetRegion_YXZ(3,:)));

        %%
        DWI_Data = double(DWI_reg.Data);

        DWI_Data_new = zeros(size(DWI_Data));
        intensities = [sum(sum(sum(DWI_reg.Data(:,:,:,1)))) sum(sum(sum(DWI_reg.Data(:,:,:,2)))) ...
            sum(sum(sum(DWI_reg.Data(:,:,:,3)))) sum(sum(sum(DWI_reg.Data(:,:,:,4))))];
        
        [intensities_sorted,I] = sort(intensities,'descend')
        for n=1:4
            DWI_Data_new(:,:,:,n) = DWI_Data(:,:,:,I(n));
        end
        DWI_Data = DWI_Data_new;

        [ysize, xsize, zsize, numbval] = size(DWI_Data);
        %bval     = DWI.bval;
        bval = [0 100 600 800];
        
        %% ADC fitting
        numvox = ysize*xsize*zsize;
        DWI_logData = log(reshape(DWI_Data,[numvox,numbval]));
        R = zeros(numvox,1);
        B = zeros(numvox,1);
        ADC = zeros(numvox,1);
        Fitting = zeros(size(DWI_logData));
        
        for ii=1:numvox
            if ~mod(ii, 500000)
                disp([num2str(ii) ' of ' num2str(numvox)]);
            end
            
            logData_ind = DWI_logData(ii,:);
            
            
            if sum(diff(logData_ind)>0)>0 || sum(isinf(logData_ind))>0
                R(ii) = 0; ADC(ii) = -1; B(ii) = 0;
                Fitting(ii,:) = NaN * ones(size(logData_ind));
            else
                [R(ii),ADC(ii),B(ii)] = regression(-bval,logData_ind);
                Fitting(ii,:) = ADC(ii)*(-bval) + B(ii);
            end
            
        end
        
        ADC_map = reshape(ADC, [ysize, xsize, zsize]);
        
        
        ADC_map(ADC_map < 0) = 0;
        ADC_map = ADC_map * 1e6;

        %% Nearest Neighbor Filling
        disp("Applying nearest neighbor filling")
        ADC_calc_fill = ADC_map;
        for z=2:size(ADC_calc_fill,3)-1
            for x=2:size(ADC_calc_fill,2)-1
                for y=2:size(ADC_calc_fill,1)-1
                    if ADC_calc_fill(y,x,z)==0
                        neighbor_list = [ADC_calc_fill(y+1,x,z) ADC_calc_fill(y,x+1,z) ADC_calc_fill(y,x,z+1) ...
                                         ADC_calc_fill(y-1,x,z) ADC_calc_fill(y,x-1,z) ADC_calc_fill(y,x,z-1) ...
                                         ADC_calc_fill(y+1,x+1,z) ADC_calc_fill(y+1,x,z+1) ADC_calc_fill(y,x+1,z+1) ...
                                         ADC_calc_fill(y-1,x-1,z) ADC_calc_fill(y-1,x,z-1) ADC_calc_fill(y,x-1,z-1) ...
                                         ADC_calc_fill(y+1,x-1,z) ADC_calc_fill(y+1,x,z-1) ADC_calc_fill(y,x+1,z-1) ...
                                         ADC_calc_fill(y-1,x+1,z) ADC_calc_fill(y-1,x,z+1) ADC_calc_fill(y,x-1,z+1) ...
                                         ADC_calc_fill(y+1,x+1,z-1) ADC_calc_fill(y+1,x-1,z+1) ADC_calc_fill(y-1,x+1,z+1) ...
                                         ADC_calc_fill(y-1,x-1,z+1) ADC_calc_fill(y-1,x+1,z-1) ADC_calc_fill(y+1,x-1,z-1) ...
                                         ADC_calc_fill(y+1,x+1,z+1) ADC_calc_fill(y-1,x-1,z-1)];
                        neighbor_list = nonzeros(neighbor_list);
                        if ~isempty(neighbor_list)
                            ADC_calc_fill(y,x,z) = mean(neighbor_list);
                        end
                    end
                end
            end
        end
        Calculated_ADC_reg = ADC_calc_fill;

        %ADC_struct.ADC = ADC;
        %ADC_struct.rMap = rMap;
        %ADC_struct.ov = ov;
        %Provided_ADC_reg = ParaMap_ADC_reg.Data;
        %Calculated_ADC_reg = ADC_map;
        
        save([datapath CaseNum '_Calculated_ADC_3DRigidReg.mat'],'Calculated_ADC_reg')
    


    end
end
