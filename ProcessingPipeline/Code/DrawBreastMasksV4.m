function [scan1,scan2,scan3,ManualParas] = DrawBreastMasksV4(scan1_0,scan2_0,scan3_0,BreastMaskFlag,pt)   

breastMask = zeros(size(scan1_0.roi));
scan1_0.maskbreast = zeros(size(scan1_0.roi));
scan2_0.maskbreast = zeros(size(scan1_0.roi));
scan3_0.maskbreast = zeros(size(scan1_0.roi));
ManualParas = struct();
breastfile0 = ['/Volumes/SanDiskSSD/ISPY2/Data/Processed/',pt,'/IntervisitRegistered/',pt,'_regData_Breast.mat'];

if BreastMaskFlag == 1
    %button = questdlg('Please select the file with previous breast masks.','Waiting...','Continue','Continue');
    %[breastfile, where] = uigetfile('*.mat');
    clear A
    %A = load([where breastfile],'scan2');
    A = load(breastfile0,'scan2');
    scan1_0.maskbreast = A.scan2.maskbreast;
    scan2_0.maskbreast = A.scan2.maskbreast;
    scan3_0.maskbreast = A.scan2.maskbreast;
elseif BreastMaskFlag == 2
    %%%%% step 1: get the breast region of pre-constrast anatomical image  
    load(['/Volumes/SanDiskSSD/ISPY2/Data/Processed/',pt,'/IntervisitRegistered/',pt,'_regData_Breast.mat'],'scan2');

    C = scan2_0.enhanced; %+scan1.enhanced+scan3.enhanced;
    
    % D = sum(C(:,:,:),3);              
    % figure; 
    % set(gcf, 'Position', get(0, 'Screensize'));
    % imagesc(D); 
    % colormap gray; 
    % axis equal; 
    % axis off;
    % title('Select breast region')
    % [ROI,~,~] = roipoly;

    ROI = scan2.maskbreast;
    box = regionprops(ROI,'BoundingBox');
    close all

    BreastRegion_YXZ = [floor(box.BoundingBox(2)) floor(box.BoundingBox(2))+ceil(box.BoundingBox(4)); ...
                        floor(box.BoundingBox(1)) floor(box.BoundingBox(1))+ceil(box.BoundingBox(3)); ...
                        1 size(scan2_0.roi,3)];
        
    ManualParas.BreastRegion_YXZ = BreastRegion_YXZ;
    
    %%%%% step 2: smooth to remove image noise 
    close all
    Anatomical_smth = zeros(size(C));
    for slice = 1:size(C,3)
        Anatomical_smth(:,:,slice) = imgaussfilt(C(:,:,slice), 5);
    end
    
    temp = zeros(size(C));
    temp(BreastRegion_YXZ(1,1):BreastRegion_YXZ(1,2), ...
         BreastRegion_YXZ(2,1):BreastRegion_YXZ(2,2), ...
         BreastRegion_YXZ(3,1):BreastRegion_YXZ(3,2)) ...
         = Anatomical_smth(BreastRegion_YXZ(1,1):BreastRegion_YXZ(1,2), ...
                           BreastRegion_YXZ(2,1):BreastRegion_YXZ(2,2), ...
                           BreastRegion_YXZ(3,1):BreastRegion_YXZ(3,2));
    Anatomical_smth = 255 * temp / prctile(temp(:),95);
    
    
    %%%%% step 3: thresholding and add tissue bounds!
    vuOnePaneViewer(Anatomical_smth)
    
    prompt = 'Breast threshold? \n';
    x = input(prompt);
    ManualParas.BreastThreshold = x;
    
    BreastMask_Int = (Anatomical_smth > x);
    BreastMask_Int(BreastRegion_YXZ(1,2),:,:) = 1;
    BreastMask_Int(BreastRegion_YXZ(1,1),:,:) = 1;
   
    %%%%% step 4: fill the holes and remove "bounds" added in step 3
    for z = 1:size(BreastMask_Int,3)
        BreastMask_Int(:,:,z) = imfill(BreastMask_Int(:,:,z), 'holes');
    end
    
    BreastMask_Int(BreastRegion_YXZ(1,1),:,:) = 0;
    BreastMask_Int(:,BreastRegion_YXZ(2,2),:) = 0;
    
    for z = 1:size(BreastMask_Int,3)
        BreastMask_Int(:,:,z) = imfill(BreastMask_Int(:,:,z), 'holes');
    end
    
    %%%%% step 5: adjust the countour
    se = strel('disk',3);
    BreastMask = zeros(size(BreastMask_Int));
    for z = 1:size(BreastMask_Int,3)
        BreastMask(:,:,z) = imerode(BreastMask_Int(:,:,z), se);
        gg = round(0.5*sum(BreastMask(:,:,z),'all'));
        BreastMask(:,:,z) = bwareaopen(BreastMask(:,:,z),gg);
    end
    
    %%%% step 6: visualization -- observation
    % for z = 1:size(BreastMask,3) 
    %     subplot(1,3,1)
    %     imagesc(C(:,:,z) .* BreastMask(:,:,z))
    %     subplot(1,3,2)
    %     imagesc(BreastMask(:,:,z))
    %     subplot(1,3,3)
    %     imagesc(Anatomical_smth(:,:,z))
    %     caxis([0 255])
    %     title([num2str(z) ' of ' num2str(size(BreastMask,3))])
    %     pause(0.1)
    % end

    close all
    
    %%%%% step 7: manual modification if necessary: mostly for the first and
    %%%%% last few slices
    % x = [];
    % while isempty(x)
    %     prompt = 'Remove any begining slices? \n Type NaN to remove none, \n otherwise type the number of slices to remove from 1 to X. \n';
    %     x = input(prompt);
    % end
    slices = find(sum(sum(scan2.roi))>0);

    x = slices(1)-1;
    if x > 0 
        BreastMask(:,:,1:x) = 0;
    end
    
    % y = [];
    % while isempty(y)
    %     prompt = 'Remove end slices? \n Type NaN to remove none, \n otherwise type the slice to remove from X to the end. \n';
    %     y = input(prompt);
    % end
    y = slices(end)+1;
    if y < size(BreastMask,3)+1
        BreastMask(:,:,y:end) = 0;
    end
    
    ManualParas.RemoveSlices = [1:x, y:size(BreastMask,3)];
        
    %%%%% step 8: smooth mask if necessary: performed to avoid manual
    %%%%% modification in step 7 introducing discontinuity between slices
    BreastMaskSmooth = double(BreastMask);
    for z = 1:size(BreastMaskSmooth,3)
        BreastMaskSmooth(:,:,z) = imgaussfilt3(BreastMaskSmooth(:,:,z), 4);
    end
    BreastMaskSmooth = double(BreastMaskSmooth>0.5);

    temp = zeros(size(BreastMaskSmooth,2),size(BreastMaskSmooth,3),size(BreastMaskSmooth,1));
    for y = 1:size(BreastMaskSmooth,1)
        temp(:,:,y) = imgaussfilt3(permute(BreastMaskSmooth(y,:,:),[2 3 1]), 4);
    end
    BreastMaskSmooth = permute(double(temp>0.5),[3 1 2]);
    BreastMask = BreastMaskSmooth;
    
    BreastMaskSmooth = imgaussfilt3(BreastMaskSmooth,4);
    BreastMaskSmooth = double(BreastMaskSmooth>0.5);
    
    
    for z = 1:size(BreastMask,3)
        gg = round(0.5*sum(BreastMask(:,:,z),'all'));
        BreastMask(:,:,z) = bwareaopen(BreastMask(:,:,z),gg);
    end

    for z = 1:size(BreastMaskSmooth,3) 
        subplot(1,3,1)
        imagesc(C(:,:,z) .* BreastMask(:,:,z))
        subplot(1,3,2)
        imagesc(BreastMask(:,:,z))
        subplot(1,3,3)
        imagesc(Anatomical_smth(:,:,z))
        title(z)
        pause(0.1)
    end
        
    close all

    %%% step 9: save mask
    fprintf('Saving breast masks\n');

    scan1_0.maskbreast = BreastMask;
    scan2_0.maskbreast = BreastMask;
    scan3_0.maskbreast = BreastMask;
else
    %For each slice draw breast mask and keep track of max and min mask
    %dimensions
%     vuThreePaneViewer(scan2.avgdce+scan1.avgdce+scan3.avgdce...
%                       +100*(scan2.roi+scan1.roi+scan3.roi));

    %%%%% step 1: get the breast region of pre-constrast anatomical image              
    C = scan2_0.enhanced; %+scan1.enhanced+scan3.enhanced;
    
    D = sum(C(:,:,:),3);              
    figure; 
    set(gcf, 'Position', get(0, 'Screensize'));
    imagesc(D); 
    colormap gray; 
    axis equal; 
    axis off;
    title('Select breast region')
    [ROI,~,~] = roipoly;
    box = regionprops(ROI,'BoundingBox');
    close all

    BreastRegion_YXZ = [floor(box.BoundingBox(2)) floor(box.BoundingBox(2))+ceil(box.BoundingBox(4)); ...
                        floor(box.BoundingBox(1)) floor(box.BoundingBox(1))+ceil(box.BoundingBox(3)); ...
                        1 size(scan2_0.roi,3)];
        
    ManualParas.BreastRegion_YXZ = BreastRegion_YXZ;
    
    %%%%% step 2: smooth to remove image noise                                
    Anatomical_smth = zeros(size(C));
    for slice = 1:size(C,3)
        Anatomical_smth(:,:,slice) = imgaussfilt(C(:,:,slice), 5);
    end
    
    temp = zeros(size(C));
    temp(BreastRegion_YXZ(1,1):BreastRegion_YXZ(1,2), ...
         BreastRegion_YXZ(2,1):BreastRegion_YXZ(2,2), ...
         BreastRegion_YXZ(3,1):BreastRegion_YXZ(3,2)) ...
         = Anatomical_smth(BreastRegion_YXZ(1,1):BreastRegion_YXZ(1,2), ...
                           BreastRegion_YXZ(2,1):BreastRegion_YXZ(2,2), ...
                           BreastRegion_YXZ(3,1):BreastRegion_YXZ(3,2));
    Anatomical_smth = 255 * temp / prctile(temp(:),95);
    
    
    %%%%% step 3: thresholding and add tissue bounds!
    vuOnePaneViewer(Anatomical_smth)
    
    prompt = 'Breast threshold? \n';
    x = input(prompt);
    ManualParas.BreastThreshold = x;
    
    BreastMask_Int = (Anatomical_smth > x);
    BreastMask_Int(BreastRegion_YXZ(1,2),:,:) = 1;
    BreastMask_Int(BreastRegion_YXZ(1,1),:,:) = 1;
    %BreastMask_Int(:,BreastRegion_YXZ(2,2),:) = 1;
    %BreastMask_Int(:,BreastRegion_YXZ(2,1),:) = 1;

    %%%%% step 3.1: identify chest wall region
    %{
    D = sum(Anatomical_smth,3);              
    figure; 
    set(gcf, 'Position', get(0, 'Screensize'));
    imagesc(D); 
    colormap gray; 
    axis equal; 
    axis off;
    title('Select chest wall region')
    [ROI,~,~] = roipoly;
    box = regionprops(ROI,'BoundingBox');
    close all

    ChestRegion_YXZ = [BreastRegion_YXZ(1,1) floor(box.BoundingBox(2))+ceil(box.BoundingBox(4)); ...
                       BreastRegion_YXZ(2,1) floor(box.BoundingBox(1))+ceil(box.BoundingBox(3)); ...
                       1 size(scan2.roi,3)];
        
    ManualParas.ChestRegion_YXZ = ChestRegion_YXZ;
    
    %%%%% step 3.2: thresholding to remove chest wall
    temp = zeros(size(Anatomical_smth));
    temp(ChestRegion_YXZ(1,1):ChestRegion_YXZ(1,2), ...
         ChestRegion_YXZ(2,1):ChestRegion_YXZ(2,2), ...
         ChestRegion_YXZ(3,1):ChestRegion_YXZ(3,2))...
         = Anatomical_smth(ChestRegion_YXZ(1,1):ChestRegion_YXZ(1,2), ...
                           ChestRegion_YXZ(2,1):ChestRegion_YXZ(2,2), ...
                           ChestRegion_YXZ(3,1):ChestRegion_YXZ(3,2));
    vuOnePaneViewer(temp)
    prompt = 'Chest wall threshold? \n';
    x = input(prompt);
    ManualParas.ChestWallThreshold = x;
    
    temp = (temp > x);
    temp(ChestRegion_YXZ(1,2),:,:) = 1;  
    temp(:,ChestRegion_YXZ(2,2),:) = 1;
    temp(:,ChestRegion_YXZ(2,1),:) = 1;

    for z = 1:size(temp,3)
        temp(:,:,z) = imfill(temp(:,:,z), 'holes');
    end
    temp(ChestRegion_YXZ(1,2),:,:) = 0;  
    temp(:,ChestRegion_YXZ(2,2),:) = 0;
    temp(:,ChestRegion_YXZ(2,1),:) = 0;
    %}
    
    %%%%% step 4: fill the holes and remove "bounds" added in step 3
    for z = 1:size(BreastMask_Int,3)
        BreastMask_Int(:,:,z) = imfill(BreastMask_Int(:,:,z), 'holes');
    end
    
    BreastMask_Int(BreastRegion_YXZ(1,1),:,:) = 0;
    BreastMask_Int(:,BreastRegion_YXZ(2,2),:) = 0;
    
% % %     BreastMask_Int(temp) = 0;
    for z = 1:size(BreastMask_Int,3)
        BreastMask_Int(:,:,z) = imfill(BreastMask_Int(:,:,z), 'holes');
    end
    
    %%%%% step 5: adjust the countour
    se = strel('disk',3);
    BreastMask = zeros(size(BreastMask_Int));
    for z = 1:size(BreastMask_Int,3)
        BreastMask(:,:,z) = imerode(BreastMask_Int(:,:,z), se);
        gg = round(0.5*sum(BreastMask(:,:,z),'all'));
        BreastMask(:,:,z) = bwareaopen(BreastMask(:,:,z),gg);
    end
    
    %%%% step 6: visualization -- observation
    for z = 1:size(BreastMask,3) 
        subplot(1,3,1)
        imagesc(C(:,:,z) .* BreastMask(:,:,z))
        subplot(1,3,2)
        imagesc(BreastMask(:,:,z))
        subplot(1,3,3)
        imagesc(Anatomical_smth(:,:,z))
        caxis([0 255])
        title([num2str(z) ' of ' num2str(size(BreastMask,3))])
        pause(0.1)
    end

    close all
    
    %%%%% step 7: manual modification if necessary: mostly for the first and
    %%%%% last few slices
    x = [];
    while isempty(x)
        prompt = 'Remove any begining slices? \n Type NaN to remove none, \n otherwise type the number of slices to remove from 1 to X. \n';
        x = input(prompt);
    end
    if x > 0
        BreastMask(:,:,1:x) = 0;
    end
    
    y = [];
    while isempty(y)
        prompt = 'Remove end slices? \n Type NaN to remove none, \n otherwise type the slice to remove from X to the end. \n';
        y = input(prompt);
    end
    if y > 0
        BreastMask(:,:,y:end) = 0;
    end
    
    ManualParas.RemoveSlices = [1:x, y:size(BreastMask,3)];
    
        
    %%%%% step 8: smooth mask if necessary: performed to avoid manual
    %%%%% modification in step 7 introducing discontinuity between slices
    BreastMaskSmooth = double(BreastMask);
    for z = 1:size(BreastMaskSmooth,3)
        BreastMaskSmooth(:,:,z) = imgaussfilt3(BreastMaskSmooth(:,:,z), 4);
    end
    BreastMaskSmooth = double(BreastMaskSmooth>0.5);

    temp = zeros(size(BreastMaskSmooth,2),size(BreastMaskSmooth,3),size(BreastMaskSmooth,1));
    for y = 1:size(BreastMaskSmooth,1)
        temp(:,:,y) = imgaussfilt3(permute(BreastMaskSmooth(y,:,:),[2 3 1]), 4);
    end
    BreastMaskSmooth = permute(double(temp>0.5),[3 1 2]);
    BreastMask = BreastMaskSmooth;
    
    BreastMaskSmooth = imgaussfilt3(BreastMaskSmooth,4);
    BreastMaskSmooth = double(BreastMaskSmooth>0.5);
    
    
    for z = 1:size(BreastMask,3)
        gg = round(0.5*sum(BreastMask(:,:,z),'all'));
        BreastMask(:,:,z) = bwareaopen(BreastMask(:,:,z),gg);
    end

    for z = 1:size(BreastMaskSmooth,3) 
        subplot(1,3,1)
        imagesc(C(:,:,z) .* BreastMask(:,:,z))
        subplot(1,3,2)
        imagesc(BreastMask(:,:,z))
        subplot(1,3,3)
        imagesc(Anatomical_smth(:,:,z))
        title(z)
        pause(0.1)
    end
        
    close all

    %%% step 9: save mask
    fprintf('Saving breast masks\n');

    scan1_0.maskbreast = BreastMask;
    scan2_0.maskbreast = BreastMask;
    scan3_0.maskbreast = BreastMask;
end

scan1 = scan1_0;
scan2 = scan2_0;
scan3 = scan3_0;
