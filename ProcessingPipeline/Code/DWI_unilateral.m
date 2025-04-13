clear
clc
startpath = pwd;

% % matpath = '/net/bananastand/data1/MDACC_breast/Processed/';
% % matpath    = '/Volumes/SanDiskSSD/ISPY2/Data/Processed/';
matpath    = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';


VariableNames = {'DWI','ParaMap_ADC','ParaMap_DWImask'}; %

Patient_list = {
'104268'...
};

%% %%% for convinience of post-processing, all formated data are moved to a
%%%%% single folder

for patientID = 1:numel(Patient_list)
    PatientName = Patient_list{patientID};
    % matpath  = [matpath,   PatientName,  '/UnregisteredData/'];
    savpath = [matpath, PatientName, '/'];
    figpath  = [matpath, PatientName, '/figures/'];
    if ~isfolder(figpath)
        mkdir(figpath)
    end
    
    for scanID = 1:2 
        CaseNum = [PatientName 'T' num2str(scanID)]

        for vidx = 1:length(VariableNames)
            load([matpath, CaseNum, '_',VariableNames{vidx},'_scandata.mat']);
        end
        % load([matpath, CaseNum, '_LocInfo.mat']);
        
        %% visual & draw
        close all
        figure(1)
        imagesc(max(DWI.Data(:,:,:,end), [], 3))
            
            BreastBound = roipoly();
            
        %%   
        DWI0 = DWI;

        YRange = find(sum(BreastBound, 2)); %13:115
        XRange = find(sum(BreastBound, 1)); %103:187
         
        pos0_Y = DWI0.pos(1,2,1) + (YRange(1) - 1) * DWI0.resolution(1);
        pos0_X = DWI0.pos(1,1,1) + (XRange(1) - 1) * DWI0.resolution(2);
        row     = length(YRange);
        columns = length(XRange);
        
        for vidx = 1:length(VariableNames)
            varname = VariableNames{vidx};
            % eval([varname, ' = ',varname,';'])
            eval([varname, '.Data = ',varname,'.Data(YRange, XRange, :, :);'])
            eval([varname, '.pos(:,2,:) = pos0_Y;'])
            eval([varname, '.pos(:,1,:) = pos0_X;'])
            eval([varname, '.row        = row;'])
            eval([varname, '.columns    = columns;'])
            eval([varname, '.TrimRange_YX = [YRange(1), XRange(1); YRange(end), XRange(end);];'])
            
            save([savpath, CaseNum, '_',varname,'_scandata.mat'],varname,'-v7.3');
        end
    end
end