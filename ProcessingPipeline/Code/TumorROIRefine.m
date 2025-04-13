function [scan1, scan2, scan3] = TumorROIRefine(scan1, scan2, scan3)
    ROI1 = scan1.roi;
    ROI2 = scan2.roi;
    ROI3 = scan3.roi;

    se1 = strel('disk',2);
%     se1 = strel('disk',1);
    %%
    ROI1_temp = imdilate(ROI1,se1);
    slices = squeeze(find(sum(sum(ROI1_temp))));
    for idx = 1:length(slices)
        i = slices(idx);
        ROI1_temp(:,:,i) = imfill(ROI1_temp(:,:,i),'holes');
    end
    ROI1_temp1 = imerode(ROI1_temp,se1);
    ROI1_temp2 = double(imgaussian(ROI1_temp1,1) >0.5);
    
    ROI2_temp = imdilate(ROI2,se1);
    slices = squeeze(find(sum(sum(ROI2_temp))));
    for idx = 1:length(slices)
        i = slices(idx);
        ROI2_temp(:,:,i) = imfill(ROI2_temp(:,:,i),'holes');
    end
    ROI2_temp1 = imerode(ROI2_temp,se1);
    ROI2_temp2 = double(imgaussian(ROI2_temp1,1) >0.5);
    
    
    ROI3_temp = imdilate(ROI3,se1);
    slices = squeeze(find(sum(sum(ROI3_temp))));
    for idx = 1:length(slices)
        i = slices(idx);
        ROI3_temp(:,:,i) = imfill(ROI3_temp(:,:,i),'holes');
    end
    ROI3_temp1 = imerode(ROI3_temp,se1);
    ROI3_temp2 = double(imgaussian(ROI3_temp1,1) >0.5);
    
    %%
    scan1.roi      = ROI1_temp2; 
    scan2.roi      = ROI2_temp2; 
    scan3.roi      = ROI3_temp2; 
    
    