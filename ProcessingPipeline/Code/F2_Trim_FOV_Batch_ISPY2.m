clear
clc
startpath = pwd;

% % matpath = '/net/bananastand/data1/MDACC_breast/Processed/';
% % matpath    = '/Volumes/SanDiskSSD/ISPY2/Data/Processed/';
matpath    = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';


VariableNames = {'DWI','DCE','ParaMap_ADC','ParaMap_DWImask'}; %

Patient_list = {
'104268'...
};

%% %%% for convinience of post-processing, all formated data are moved to a
%%%%% single folder

for patientID = 1:numel(Patient_list)
    PatientName = Patient_list{patientID};
    % matpath  = [matpath,   PatientName,  '/UnregisteredData/'];
    
    figpath  = [matpath, PatientName, '/figures/'];
    if ~isfolder(figpath)
        mkdir(figpath)
    end
    
    for scanID = 0 %:2 
        
        CaseNum = [PatientName 'T' num2str(scanID)]

        for vidx = 1:length(VariableNames)
            load([matpath, CaseNum, '_',VariableNames{vidx},'_align.mat']);
        end
        load([matpath, CaseNum, '_LocInfo.mat']);
        
        %% visual & draw
        close all
        figure(1)
        imagesc(max(DCE_align.Data(:,:,:,end) .* (1 + ParaMap_DWImask_align.Data), [], 3))
            
            BreastBound = roipoly();
            
        %%     
        YRange = find(sum(BreastBound, 2));
        XRange = find(sum(BreastBound, 1));
         
        pos0_Y = DCE_align.pos(1,1,1) + (YRange(1) - 1) * DCE_align.resolution(1);
        pos0_X = DCE_align.pos(1,2,1) + (XRange(1) - 1) * DCE_align.resolution(2);
        row     = length(YRange);
        columns = length(XRange);
        
        for vidx = 1:length(VariableNames)
            varname = VariableNames{vidx};
            eval([varname, '_trim = ',varname,'_align;'])
            eval([varname, '_trim.Data = ',varname,'_align.Data(YRange, XRange, :, :);'])
            eval([varname, '_trim.pos(:,1,:) = pos0_Y;'])
            eval([varname, '_trim.pos(:,2,:) = pos0_X;'])
            eval([varname, '_trim.row        = row;'])
            eval([varname, '_trim.columns    = columns;'])
            eval([varname, '_trim.TrimRange_YX = [YRange(1), XRange(1); YRange(end), XRange(end);];'])
            
            save([matpath, CaseNum, '_',varname,'_trim.mat'],[varname,'_trim'],'-v7.3');
        end
         
        
    end
end
