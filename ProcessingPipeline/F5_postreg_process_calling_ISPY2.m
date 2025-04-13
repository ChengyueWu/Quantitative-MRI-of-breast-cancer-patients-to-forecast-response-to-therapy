clear 
close all
clc;

%Set flag to 1 if breast masks have already been drawn and you will be
%prompted to select the 'RAW' file containing the breast masks before
%trimming

Scriptpath = [pwd,'/'];
filepath = '/Volumes/SanDiskSSD/ISPY2/Data/Processed/';

addpath(Scriptpath)
addpath(filepath) 

%{
'698094'...
'831848'...
'414844'...
'991818'...
'896852'...

%}

Patient_list = {
'414844'...
};

BreastMaskFlag    = 1; % 1 = already saved, 2 = need to draw
filenameSubfix    = '_regData_Breast';
savenameSubfix    = '_regData';

targetimageforreg = 2; %Default is 2

%%                    
for patientID = 1:numel(Patient_list)
    tic
    subname = Patient_list{patientID};
    subname

    figpath  = [filepath, subname, '/figures/'];
    if ~isfolder(figpath)
        mkdir(figpath)
    end
%
%% Make structures of the registered data and return image sizes
[imagesizes,scan1,scan2,scan3] = GenerateScanStructs_postReg_ISPY2(filepath,subname);

%% Draw breast masks -- Manual input required!!
fprintf('Starting breast mask drawer. \n');
[scan1,scan2,scan3, MPs_BreastMask] = DrawBreastMasksV4(scan1,scan2,scan3,BreastMaskFlag,subname);
fprintf('Breast masks generated. \n');

    %%%% plot (optional)
    figure
    col=[0.3 0.3 0.3];
    hiso1 = patch(isosurface(scan1.maskbreast,0),'FaceColor',col,'EdgeColor','none');
    hiso2 = patch(isocaps(scan1.maskbreast,0),   'FaceColor',col,'EdgeColor','none');
    alpha(0.3);
    
    grid on
    lighting phong;
    camlight;
    
    %% check if further manual modification is needed ---
    %{
    % load...
    filepath = '/Volumes/ChengyueWu2/DataBackUp/MDACC_ClinicalPipeline/MDACC_patient_processing/Processed/';
    for z = 1:size(scan2.maskbreast,3) 
        subplot(1,3,1)
        imagesc(scan2.avgdce(:,:,z) .* scan2.maskbreast(:,:,z))
        subplot(1,3,2)
        imagesc(scan2.maskbreast(:,:,z))
        subplot(1,3,3)
        imagesc(scan2.avgdce(:,:,z))
        title(z)
        pause
    end
%     imagesc(max(permute(scan2.maskbreast,[1 3 2]),[],3))
    
    
    %% "-" ---
        close all
        MaskMultiply = ones(size(scan2.maskbreast));
        
        figure
        ManualAdj_slices = [27:33, 113:121];
        for z = ManualAdj_slices
            imagesc(scan2.avgdce(:,:,z) .* scan2.maskbreast(:,:,z))
            title(z)
            MaskMultiply(:,:,z) = roipoly();
        end
        MPs_BreastMask.MaskMultiply = MaskMultiply;
        Mask = double(imgaussfilt3(MaskMultiply .* scan2.maskbreast,4) > 0.5);
        
        scan1.maskbreast = Mask;
        scan2.maskbreast = Mask;
        scan3.maskbreast = Mask;
        
        
% %         LocInfo = load([filepath, subname ,'/UnregisteredData/',subname '_v1_LocInfo.mat']); %
% %         dx = mean(abs(diff(LocInfo.Xlist_dce)));
% %         dy = mean(abs(diff(LocInfo.Ylist_dce)));
% %         dz = mean(abs(diff(LocInfo.Zlist_dce)));
% %         dims = [dx, dy, dz]
% %         
% %         savename = [subname '_regData_Breast']
% %         filepath = '/Volumes/ChengyueWu2/DataBackUp/MDACC_ClinicalPipeline/MDACC_patient_processing/Processed/';
% %         save([filepath subname, '/IntervisitRegistered/', savename '.mat'],...
% %         'scan1','scan2','scan3','dims','dx','dy','dz','subname','MPs_BreastMask','-v7.3');

    %%  "+" ---
        close all
        MaskAdd = zeros(size(scan2.maskbreast));
        
        figure
        ManualAdj_slices = [117:121];
        for z = ManualAdj_slices
            imagesc(scan2.avgdce(:,:,z) .* (0.3+scan2.maskbreast(:,:,z)))
            title(z)
            MaskAdd(:,:,z) = roipoly();
        end
        MPs_BreastMask.MaskAdd = MaskAdd;
        Mask = double(imgaussfilt3(double(logical(MaskAdd + scan2.maskbreast)),4) > 0.5);
        
        scan1.maskbreast = Mask;
        scan2.maskbreast = Mask;
        scan3.maskbreast = Mask;
        
        
        
        LocInfo = load([filepath, subname ,'/UnregisteredData/',subname '_v1_LocInfo.mat']); %
        dx = mean(abs(diff(LocInfo.Xlist_dce)));
        dy = mean(abs(diff(LocInfo.Ylist_dce)));
        dz = mean(abs(diff(LocInfo.Zlist_dce)));
        dims = [dx, dy, dz]

        savename = [subname '_regData_Breast']
        filepath = '/Volumes/ChengyueWu2/DataBackUp/MDACC_ClinicalPipeline/MDACC_patient_processing/Processed/';
        save([filepath subname, '/IntervisitRegistered/', savename '.mat'],...
        'scan1','scan2','scan3','dims','dx','dy','dz','subname','MPs_BreastMask','-v7.3');
    
        %}
        
%% save breast segmentation 
%
LocInfo = load([filepath, subname ,'/UnregisteredData/',subname 'T0_LocInfo.mat']); %
dx = mean(abs(diff(LocInfo.Xlist_dce)));
dy = mean(abs(diff(LocInfo.Ylist_dce)));
dz = mean(abs(diff(LocInfo.Zlist_dce)));
dims = [dx, dy, dz]

savename = [subname '_regData_Breast']
save([filepath subname, '/IntervisitRegistered/', savename '.mat'],...
     'scan1','scan2','scan3','dims','dx','dy','dz','subname','MPs_BreastMask','-v7.3');

cd(Scriptpath)

%}

%% tissue segments
%{
%Set flag to 1 if tissue masks have already been drawn and you will be prompted to select the 'RAW' file containing the breast masks before trimming
TissueMaskFlag    = 0;

filename = [subname filenameSubfix]; %_MA
load([filepath subname, '/IntervisitRegistered/', filename '.mat']);

%%%%%% Generate masks for adipose and fibro tissues    
fprintf('Segmenting adipose and fibroglandular tissues. \n');
[scan1,scan2,scan3] = BreastTissuesSegmentation(scan1,scan2,scan3,TissueMaskFlag); %MPs_TissueMask
fprintf('Tissue maps generated. \n');

%%%%%% Saving data prior to trimming and smoothing
fprintf('Saving data prior to trimming and smoothing. \n');
%Saving all the slices
savename = [subname '_regData']; 
savpath = [filepath '/' subname, '/IntervisitRegistered/'];
save([savpath savename '.mat'],...
     'scan1','scan2','scan3','dims','dx','dy','dz','subname','MPs_BreastMask' , '-v7.3');
fprintf(['File saved as ' savename '.mat \n\n\n'])

%% calculations
%%% Calculate NTCs
scan1.roi = double(scan1.roi);
scan2.roi = double(scan2.roi);
scan3.roi = double(scan3.roi);

fprintf('Calculating tumor cells and data bounds. \n');
[scan1,scan2,scan3,cellradius,cellVolume,voxVolume,CarryCap,ADCw,ADCmin] = CalculateNTCs(scan1,scan2,scan3,dims);       

%%% Calculate sizes of tumors
fprintf('Calculating volume, total cells, and longest axis. \n');
[scan1,scan2,scan3] = SizingTumors(scan1,scan2,scan3,voxVolume,dims);

%%% Determine smallest domain
fprintf('Determining smallest modeling domain. \n');
[X,Y] = TrimmingWindow(scan2);

%%% Saving data prior to trimming
fprintf('Saving data prior to trimming. \n');

%Saving all the slices
savename = [subname '_regData']; 
savpath  = [filepath '/' subname, '/IntervisitRegistered/'];
save([savpath savename '.mat'],...
     'scan1','scan2','scan3','ADCw','ADCmin','dims','dx','dy','dz', ...
     'cellradius','cellVolume','voxVolume','CarryCap','subname','X','Y', '-v7.3');

fprintf(['File saved as ' savename '.mat \n\n\n'])

    %% save to new subfolders -- in preparation for modeling
    %
    %Defining all the data for the window view for modeling
    fprintf('Defining MATLAB files on the smaller window. \n');

    whichslices = find(sum(sum(scan1.maskbreast,1),2)>0);
    startsl = whichslices(1);
    endsl   = whichslices(end);

    %Breast Masks
    BreastMask = scan1.maskbreast(Y(:,1),X(1,:),startsl:endsl);

    %Tumor cell numbers
    NTC1 = scan1.NTC(Y(:,1),X(1,:),startsl:endsl);
    NTC2 = scan2.NTC(Y(:,1),X(1,:),startsl:endsl);
    NTC3 = scan3.NTC(Y(:,1),X(1,:),startsl:endsl);

    %ADC values
    ADC1 = scan1.adc(Y(:,1),X(1,:),startsl:endsl);
    ADC2 = scan2.adc(Y(:,1),X(1,:),startsl:endsl);
    ADC3 = scan3.adc(Y(:,1),X(1,:),startsl:endsl);

    %Anatomical
    Anatomical1 = scan1.avgdce(Y(:,1),X(1,:),startsl:endsl);
    Anatomical2 = scan2.avgdce(Y(:,1),X(1,:),startsl:endsl);
    Anatomical3 = scan3.avgdce(Y(:,1),X(1,:),startsl:endsl);

    %Anatomical with adaptive equalization
    AnatomicalEnh1 = scan1.enhanced(Y(:,1),X(1,:),startsl:endsl);
    AnatomicalEnh2 = scan2.enhanced(Y(:,1),X(1,:),startsl:endsl);
    AnatomicalEnh3 = scan3.enhanced(Y(:,1),X(1,:),startsl:endsl);

    %Breast tissue properties
    Tissues1 = 2*scan1.maskfibro(Y(:,1),X(1,:),startsl:endsl) + ...
                 scan1.maskadipo(Y(:,1),X(1,:),startsl:endsl);
    Tissues2 = 2*scan2.maskfibro(Y(:,1),X(1,:),startsl:endsl) + ...
                 scan2.maskadipo(Y(:,1),X(1,:),startsl:endsl);
    Tissues3 = 2*scan3.maskfibro(Y(:,1),X(1,:),startsl:endsl) + ...
                 scan3.maskadipo(Y(:,1),X(1,:),startsl:endsl);
    Tissues1(Tissues1>2) = 2;
    Tissues2(Tissues2>2) = 2;
    Tissues3(Tissues3>2) = 2;

    %ROI
    ROI1 = scan1.roi(Y(:,1),X(1,:),startsl:endsl);
    ROI2 = scan2.roi(Y(:,1),X(1,:),startsl:endsl);
    ROI3 = scan3.roi(Y(:,1),X(1,:),startsl:endsl);

    cd([filepath '/' subname]);
    filename = 'IntervisitRegistered';
    cd(filename) 

    filename = 'ModelMatFiles';
    if isfolder(filename) == 0
        mkdir(filename)
    end
    cd(filename) 

    %Saving all the files for modeling
    %Save native X and Y
    matname = ['NativeX_' subname '.mat']; save(matname,'X');
    matname = ['NativeY_' subname '.mat']; save(matname,'Y');

    %Save breast masks
    matname = ['BreastMask_' subname '.mat']; save(matname,'BreastMask');

    %Save Scan tumor cell numbers
    matname = ['NTC1_' subname '.mat']; save(matname,'NTC1');
    matname = ['NTC2_' subname '.mat']; save(matname,'NTC2');
    matname = ['NTC3_' subname '.mat']; save(matname,'NTC3');

    %ADC values
    matname = ['ADC1_' subname '.mat']; save(matname,'ADC1');
    matname = ['ADC2_' subname '.mat']; save(matname,'ADC2');
    matname = ['ADC3_' subname '.mat']; save(matname,'ADC3');

    %Save breast tissue info
    matname = ['Tissues1_' subname '.mat']; save(matname,'Tissues1');
    matname = ['Tissues2_' subname '.mat']; save(matname,'Tissues2');
    matname = ['Tissues3_' subname '.mat']; save(matname,'Tissues3');

    %Save anatomical DCE info
    matname = ['Anatomical1_' subname '.mat']; save(matname,'Anatomical1');
    matname = ['Anatomical2_' subname '.mat']; save(matname,'Anatomical2');
    matname = ['Anatomical3_' subname '.mat']; save(matname,'Anatomical3');

    %Save equalized anatomical DCE info
    matname = ['AnatomicalEnh1_' subname '.mat']; save(matname,'AnatomicalEnh1');
    matname = ['AnatomicalEnh2_' subname '.mat']; save(matname,'AnatomicalEnh2');
    matname = ['AnatomicalEnh3_' subname '.mat']; save(matname,'AnatomicalEnh3');
    
    matname = ['ROI1_' subname '.mat']; save(matname,'ROI1');
    matname = ['ROI2_' subname '.mat']; save(matname,'ROI2');
    matname = ['ROI3_' subname '.mat']; save(matname,'ROI3');
    %}
    
    %%
    cd(Scriptpath)
    fprintf('Processing complete. \n');
% % % %     beep
%}
    beep
    toc
end

cd(Scriptpath)

