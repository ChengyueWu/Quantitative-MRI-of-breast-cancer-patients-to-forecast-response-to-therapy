function [scan1,scan2,scan3] = BreastTissuesSegmentation(scan1,scan2,scan3,BreastMaskFlag)    

if BreastMaskFlag == 1
    button = questdlg('Please select the file with previous breast masks.','Waiting...','Continue','Continue');
    [breastfile, where] = uigetfile('*.mat');
    clear A
    A = load([where breastfile],'scan1');
    scan1.maskfibro = A.scan1.maskfibro;
    scan1.maskadipo = A.scan1.maskadipo;
    A = load([where breastfile],'scan2');
    scan2.maskfibro = A.scan2.maskfibro;
    scan2.maskadipo = A.scan2.maskadipo;
    A = load([where breastfile],'scan3');
    scan3.maskfibro = A.scan3.maskfibro;
    scan3.maskadipo = A.scan3.maskadipo;
else

%GENERATE MASKS FOR ADIPOSE AND FIBROGLANDULAR TISSUES
for jj = 1:3
    
    eval(['anatomical = scan' num2str(jj) '.enhanced;']); 
    eval(['BreastMask = scan' num2str(jj) '.maskbreast;']);
    eval(['removetum = imcomplement(scan' num2str(jj) '.roi);']); 

    img = BreastMask.*anatomical.*removetum;
    [ysize, xsize, ~] = size(BreastMask);

    k = 5;
    I_enh = ImageEnhancement(img, k);

    Max_perslice = max(max(I_enh));
    I_enh = 255 * I_enh ./ repmat(Max_perslice, [ysize, xsize, 1]);

    I_enh(isnan(I_enh)) = 0;


    IM_fibro = I_enh .* (I_enh > 155);
    IM_adipo = I_enh .* (I_enh <= 155); 

    Mask_fibro = bwareaopen(logical(IM_fibro), 10);
    Mask_adipo = logical(IM_adipo);

%For the UT data set the following removes too much of the image
%     anatomical_adipo = anatomical .* Mask_adipo;
%     Mask_list = find(Mask_adipo);
% 
%     [AdipoIDX, C] = kmeans(anatomical_adipo(Mask_list),2);
% 
%     if C(1,end)>C(2,end)
%        AdipoIDX = AdipoIDX - 1;
%     else
%        AdipoIDX = 2-AdipoIDX;
%     end
% 
%     [X, Y, Z] = size(anatomical_adipo);
%     IM_adipo = zeros(X, Y, Z);
%     IM_adipo(Mask_list) = anatomical_adipo(Mask_list) .* AdipoIDX;
%     Mask_adipo = double(logical(IM_adipo));

    eval(['scan' num2str(jj) '.maskfibro = double(Mask_fibro);']); 
    eval(['scan' num2str(jj) '.maskadipo = double(Mask_adipo);']); 

end
    
end    
