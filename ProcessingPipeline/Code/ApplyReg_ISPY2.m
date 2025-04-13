function ApplyReg_ISPY2(subname,matpath, Visits,Visit_num_f,weight)
%Apply registration to all the data maps
Visits = Visits - 1;
Visit_num_f = Visit_num_f - 1;

cd(matpath)

for Visit_num_m = Visits  
    outpath = [matpath subname 'T' num2str(Visit_num_m) 'regtoT' num2str(Visit_num_f)' '/'];
    if Visit_num_m == Visit_num_f
        filename = [subname 'T',num2str(Visit_num_f),'_FixedImagesForReg.mat'];
    else
        filename = [subname 'T',num2str(Visit_num_m),'_MovingToBeRegTo_T',num2str(Visit_num_f) '.mat'];
    end
    clear roi_orig roi
    load(filename);
    roi = roi_orig;
    
    %removing padding
    avgdce_temp     = zeros(size(roi,1),size(roi,2),size(roi,3)+2);
    anatomical_temp = zeros(size(roi,1),size(roi,2),size(roi,3)+2);
    roi_temp        = zeros(size(roi,1),size(roi,2),size(roi,3)+2);
    adc_temp        = zeros(size(roi,1),size(roi,2),size(roi,3)+2);

    avgdce_temp(:,:,2:end-1)     = avgdce;
    anatomical_temp(:,:,2:end-1) = anatomical;
    roi_temp(:,:,2:end-1)        = roi;
    adc_temp(:,:,2:end-1)        = adc;


    if Visit_num_m ~= Visit_num_f
        filename = [subname,'T',num2str(Visit_num_f),'T',num2str(Visit_num_m)];
        load([outpath filename '.mat'],'DeformFields')
        
        %% step1: rigid
        df_idx = 1;
        D   = DeformFields{df_idx}.D;
        Ref = DeformFields{df_idx}.Ref;
        avgdce_temp_df     = imwarp(avgdce_temp,     Ref, D, 'OutputView',Ref);
        anatomical_temp_df = imwarp(anatomical_temp, Ref, D, 'OutputView',Ref);
        roi_temp_df        = imwarp(roi_temp,        Ref, D, 'OutputView',Ref);
        adc_temp_df        = imwarp(adc_temp,        Ref, D, 'OutputView',Ref);
        
        %% step2: non-rigid
        for df_idx = 2:length(DeformFields)
            D = zeros(size(DeformFields{df_idx}.data,1),size(DeformFields{df_idx}.data,2),size(DeformFields{df_idx}.data,3),3);
            D(:,:,:,1) = DeformFields{df_idx}.datay;
            D(:,:,:,2) = DeformFields{df_idx}.datax;
            D(:,:,:,3) = DeformFields{df_idx}.dataz;

            D(abs(D) == Inf) = NaN;
            D(isnan(D)) = 0;

            avgdce_temp_df     = imwarp(avgdce_temp_df, D);
            anatomical_temp_df = imwarp(anatomical_temp_df, D);
            roi_temp_df        = imwarp(roi_temp_df, D);              
            adc_temp_df        = imwarp(adc_temp_df, D);
        end
        
        roi_temp_df = (roi_temp_df >= 0.5);
        note = ['Deformed images using Rigid + Bspline with penalty, registered to visit ' num2str(Visit_num_f) '. No bounds applied after registration. Weight used ' num2str(weight*100)];

        %removing padding
        avgdce     = avgdce_temp_df(:,:,2:end-1);
        anatomical = anatomical_temp_df(:,:,2:end-1);
        roi        = roi_temp_df(:,:,2:end-1);
        adc        = adc_temp_df(:,:,2:end-1);

        %%% CHANGE HERE %%%
        savename = [subname 'T',num2str(Visit_num_m) '_Reg.mat'];
        save([matpath savename],'avgdce','anatomical','roi','adc','note', '-v7.3')
    
    else
        roi_temp = (roi_temp >= 0.5);
        note = 'Target images';
        
        %Removing padding
        avgdce     = avgdce_temp(:,:,2:end-1);
        anatomical = anatomical_temp(:,:,2:end-1);
        roi        = roi_temp(:,:,2:end-1);
        adc        = adc_temp(:,:,2:end-1);

        %%% CHANGE HERE %%%
        %Alternate file save structure
        savename = [subname 'T' num2str(Visit_num_f) '_Reg.mat'];
        save([matpath savename],'avgdce','anatomical','roi','adc','note', '-v7.3')
    end
end
