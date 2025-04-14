%% example 3: inter-visit registration

clc
clear 
close all

Scriptpath = pwd;
Elastixpath = '/Users/chengyuewu/Desktop/study_work/Lab/MDACC_ClinicalPipeline/MDACC_patient_processing/Codes/Clinical-Breast-Cancer-Model-Development-master/MDACC_DataProcessing_Reports/Scripts/elastix_macosx64/';

matpath  = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';

Patient_list = {
'104268'...
};

weights = ones(length(Patient_list),1).*0.5; %0.5 can be from 0.1 to 10

%%            
for pidx = 1:numel(Patient_list)
    subname = Patient_list{pidx};
    weight = weights(pidx);
    
    figpath  = [matpath, subname, '/figures/'];
    if ~isfolder(figpath)
        mkdir(figpath)
    end
    
    cd(Scriptpath)
    inpath  = [matpath, subname, '/'];
    outpath = [matpath, subname, '/'];
    
    %% Defining list of visits for registration based on available slices
    Visits      = [1, 2, 3];
    Visit_num_f = 2; 
    target_list = [1, 3]; 
    
        %%%% manual pre-alignment (optional)
        %%%% if locations of breast in different visit vary too much,
        %%%% manual selection of interested regions would be helpful -- 
        %%%% YRange, XRange, ZRange: 3 * L matrix
    
        %% Load in data from pipeline files
        fprintf('Generating registration files. \n');
        SetupRegFiles_ISPY2(subname,inpath,outpath,Visits,Visit_num_f,weight);

        %% run registration
        CallElastixNonRigidReg_ISPY2(subname,matpath,Elastixpath,target_list,Visit_num_f)
        fprintf('Registration complete. \n');
        
        %% Apply registration to all the data maps
        fprintf('Applying registration. \n');
        ApplyReg_ISPY2(subname,outpath,Visits,Visit_num_f,weight)

            %%
            cd(Scriptpath)
end