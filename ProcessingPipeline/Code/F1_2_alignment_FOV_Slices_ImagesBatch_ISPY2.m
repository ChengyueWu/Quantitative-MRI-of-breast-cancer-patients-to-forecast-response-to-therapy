%% example 1: Intra-visit alignment - FOV alignment 


clear
clc
close all

startpath = pwd;
matpath    = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';

PatientList = {
'104268'...
};
FlippedList = [0];

%%
for pidx = 1:length(PatientList)
    flipped = FlippedList(pidx);
    PatientID = PatientList{pidx};

    PatientName = PatientList{pidx}

    datapath  =  [matpath, PatientID, '/']; 
    
    figpath  = [matpath, PatientName, '/figures/'];
    if ~isfolder(figpath)
        mkdir(figpath)
    end
    
    for v = 0 %0:2
        
        CaseNum   = [PatientID, 'T', num2str(v)]
   
        %%
        if isfile([datapath, CaseNum,'_DCE_scandata.mat'])
            DCE_flag = 1;
            load([datapath, CaseNum,'_DCE_scandata.mat']);
        else 
            DCE_flag = 0;
        end
        if isfile([datapath, CaseNum,'_DWI_scandata.mat'])
            DWI_flag = 1;
            load([datapath, CaseNum,'_DWI_scandata.mat']);
        else
            DWI_flag = 0; 
        end
        if isfile([datapath, CaseNum,'_ParaMap_ADC_scandata.mat'])
            ParaADC_flag = 1;
            load([datapath, CaseNum,'_ParaMap_ADC_scandata.mat']);
        else
            ParaADC_flag = 0;
        end

        if isfile([datapath, CaseNum,'_ParaMap_DCEmask_scandata.mat'])
            ParaDCEMask_flag = 1;
            load([datapath, CaseNum,'_ParaMap_DCEmask_scandata.mat']);
            ParaMap_DCEmask = DCE; ParaMap_DCEmask.Data = double(TumorROI); 
        else
            ParaDCEMask_flag = 0; 
        end
        if isfile([datapath, CaseNum,'_ParaMap_DWImask_whole.mat'])
            ParaDWIMask_flag = 1;
            load([datapath, CaseNum,'_ParaMap_DWImask_whole.mat']);
        else
            ParaDWIMask_flag = 0;
        end


        %% Extract and arrange location info from datasets
        [Xlist_dce, Ylist_dce, Zlist_dce, EdgeXY_dce, ~, ~] = LocationExtract_ISPY2(DCE,flipped); 
        [Xlist_dwi, Ylist_dwi, Zlist_dwi, EdgeXY_dwi, ~, ~] = LocationExtract_ISPY2(DWI,flipped); 
% %             [Xlist_t1m, Ylist_t1m, Zlist_t1m, EdgeXY_t1m, ~, ~] = LocationExtract(T1Mapping); 
        
        %%% CHANGE HERE %%%
        Zlist_dce = sort(Zlist_dce);
        Zlist_dwi = sort(Zlist_dwi);
        %%%%%%%%%%%%%%%%%%%

        %%%%%% info for aligning DWI to DCE
        [~,xidx_DWI2Dynamic_end] = min(abs(Xlist_dce - Xlist_dwi(end)));
        [~,xidx_DWI2Dynamic_1]   = min(abs(Xlist_dce - Xlist_dwi(1)));
        [~,yidx_DWI2Dynamic_end] = min(abs(Ylist_dce - Ylist_dwi(end)));
        [~,yidx_DWI2Dynamic_1]   = min(abs(Ylist_dce - Ylist_dwi(1)));
        [~,zidx_DWI2Dynamic_end] = min(abs(Zlist_dce - Zlist_dwi(end)));
        [~,zidx_DWI2Dynamic_1]   = min(abs(Zlist_dce - Zlist_dwi(1)));

        [~,xidx_Dynamic2DWI_end] = min(abs(Xlist_dwi - Xlist_dce(end)));
        [~,xidx_Dynamic2DWI_1]   = min(abs(Xlist_dwi - Xlist_dce(1)));
        [~,yidx_Dynamic2DWI_end] = min(abs(Ylist_dwi - Ylist_dce(end)));
        [~,yidx_Dynamic2DWI_1]   = min(abs(Ylist_dwi - Ylist_dce(1)));
        [~,zidx_Dynamic2DWI_end] = min(abs(Zlist_dwi - Zlist_dce(end)));
        [~,zidx_Dynamic2DWI_1]   = min(abs(Zlist_dwi - Zlist_dce(1)));
        
%             %%% CHANGE HERE %%%
%             if zidx_Dynamic2DWI_end < zidx_Dynamic2DWI_1
%                 zidx_Dynamic2DWI_end_temp = zidx_Dynamic2DWI_end;
%                 zidx_Dynamic2DWI_1_temp = zidx_Dynamic2DWI_1;
%                 zidx_Dynamic2DWI_end = zidx_Dynamic2DWI_1_temp;
%                 zidx_Dynamic2DWI_1 = zidx_Dynamic2DWI_end_temp;
%             end
%             %%%%%%%%%%%%%%%%%%%
        
        Zlist_dwi_org = Zlist_dwi(zidx_Dynamic2DWI_1:zidx_Dynamic2DWI_end);
        Zlist_dwi_dyn = Zlist_dce(zidx_DWI2Dynamic_1:zidx_DWI2Dynamic_end);
       
        %%%%%%% central tumor slice %%%%%%%%%
        tumorslices = find(sum(sum(ParaMap_DWImask_whole.Data)));
        zloc = Zlist_dwi(round(median(tumorslices))); %round(length(Zlist_dwi)/2)
        [~, slc_dwi] = min(abs(Zlist_dwi - zloc));
        [~, slc_dce] = min(abs(Zlist_dce - zloc));

        xsize = xidx_DWI2Dynamic_end - xidx_DWI2Dynamic_1 + 1;
        ysize = yidx_DWI2Dynamic_end - yidx_DWI2Dynamic_1 + 1;
        zsize = length(Zlist_dce);

        savename = [CaseNum, '_LocInfo.mat'];
        save([datapath, savename], ...
             'Xlist_dce', 'Ylist_dce', 'Zlist_dce', 'EdgeXY_dce',...
             'Xlist_dwi', 'Ylist_dwi', 'Zlist_dwi', 'EdgeXY_dwi',...
             'xidx_DWI2Dynamic_end','xidx_DWI2Dynamic_1','yidx_DWI2Dynamic_end','yidx_DWI2Dynamic_1','zidx_DWI2Dynamic_end','zidx_DWI2Dynamic_1',...
             'xidx_Dynamic2DWI_end','xidx_Dynamic2DWI_1','yidx_Dynamic2DWI_end','yidx_Dynamic2DWI_1','zidx_Dynamic2DWI_end','zidx_Dynamic2DWI_1',...
             'Zlist_dwi_org','Zlist_dwi_dyn',...
             'zloc','slc_dwi','slc_dce',... 
             'xsize','ysize','zsize')

            %%% plots (optional): FOV & Slice ovelapping
            %
            Myplot_FOVs_Overlapping  (EdgeXY_dce, EdgeXY_dwi, NaN, NaN, figpath, CaseNum) %%EdgeXY_t1m,
            Myplot_Slices_Overlapping(Zlist_dce,  Zlist_dwi,  NaN, NaN, figpath, CaseNum) %%Zlist_t1m,


            DataSets = {DCE.Data, DWI.Data};    %% t1m_breast, 
            Names = {'DCE',   'DWI'};                 %% 'T1Mapping', 
            Slices= [slc_dce, slc_dwi];               %% slc_t1m,     

            Myplot_CentralSlices_Visual (CaseNum, figpath, DataSets, Names, Slices)
            pause
            %}
            
        %% Alignment 1: Transverse 
        %%%%%%% FOV alignment %%%%%%%%%
        %%%% step 1: trim DCE images to breast FOV (overlapping with diff)
            dce_breast = DCE.Data(yidx_DWI2Dynamic_1:yidx_DWI2Dynamic_end, ...
                                  xidx_DWI2Dynamic_1:xidx_DWI2Dynamic_end, :, :);
        DCEmask_breast = ParaMap_DCEmask.Data(yidx_DWI2Dynamic_1:yidx_DWI2Dynamic_end, ...
                                              xidx_DWI2Dynamic_1:xidx_DWI2Dynamic_end, :, :);

        %%%% step 2: align DW images and maps to the resolution of dynamics
        dwi_breast     = zeros(ysize, xsize, size(DWI.Data,3), size(DWI.Data,4));
        ADC_breast     = zeros(ysize, xsize, size(DWI.Data,3));
        DWImask_breast = zeros(ysize, xsize, size(DWI.Data,3));
        for s = 1:size(DWI.Data,3)
            for f = 1:size(DWI.Data,4)
                dwi_breast(:,:,s,f) = imresize(DWI.Data(yidx_Dynamic2DWI_1:yidx_Dynamic2DWI_end,...
                                                        xidx_Dynamic2DWI_1:xidx_Dynamic2DWI_end,s,f), ...
                                               [ysize, xsize]);
            end

                ADC_breast(:,:,s) = imresize(ParaMap_ADC.Data(yidx_Dynamic2DWI_1:yidx_Dynamic2DWI_end,...
                                                              xidx_Dynamic2DWI_1:xidx_Dynamic2DWI_end,s), ...
                                             [ysize, xsize]);
                DWImask_breast(:,:,s) = imresize(ParaMap_DWImask_whole.Data(yidx_Dynamic2DWI_1:yidx_Dynamic2DWI_end,...
                                                                            xidx_Dynamic2DWI_1:xidx_Dynamic2DWI_end,s), ...
                                             [ysize, xsize]);
        end
        
        DWImask_breast = double(DWImask_breast > 0.5);
        
            %%%% plots (optional): central slices visualization & tranverse aligment comparison
            %
            DataSets = {dce_breast, dwi_breast, ADC_breast, DCEmask_breast, DWImask_breast};    %% t1m_breast, 
            Names = {'DCE',   'DWI',    'ADC',   'DCEmask',   'DWImask'};                 %% 'T1Mapping', 
            Slices= [slc_dce, slc_dwi,  slc_dwi, slc_dce,     slc_dwi];               %% slc_t1m,     

            Myplot_CentralSlices_Visual (CaseNum, figpath, DataSets, Names, Slices)
            Myplot_CentralSlices_Compare(CaseNum, figpath, DataSets, Names, Slices)

            pause
                
            %}

        %% Alignment 2: Interpolation through slices
        %%%% step 1: Use DCE as reference grid
        NumPix     = ysize*xsize;
        CentralPix = round(NumPix / 2);
        
        %%%% step 2: for DWI images and maps --
            %%%%% maps
        disp('Align result of ADC ...')
      
        ADC_breast_org_list = reshape(ADC_breast(:,:,zidx_Dynamic2DWI_1:zidx_Dynamic2DWI_end),...
                                      [NumPix, zidx_Dynamic2DWI_end - zidx_Dynamic2DWI_1 + 1]);
        ADC_breast_dyn = zeros(ysize*xsize, zsize);
        for p = 1:NumPix
            value_org = ADC_breast_org_list(p,:)';
            value_dyn = interp1(Zlist_dwi_org, value_org, Zlist_dwi_dyn);
            ADC_breast_dyn(p, zidx_DWI2Dynamic_1:zidx_DWI2Dynamic_end) = value_dyn';
        end
            %%%%%%%% plot (optional)
            %{
            value_org = ADC_breast_org_list(CentralPix,:)';
            value_dyn = ADC_breast_dyn(CentralPix,max([zidx_DWI2Dynamic_1-5, 1]):min([zidx_DWI2Dynamic_end+5, zsize]))';
            zlist_org = Zlist_dwi_org;
            zlist_dyn = Zlist_dce(max([zidx_DWI2Dynamic_1-5, 1]):min([zidx_DWI2Dynamic_end+5, zsize]));
            savename = [CaseNum, '_Slicealignment_t1m-ADC_central'];
            Myplot_SliceResample_Compare(value_org, value_dyn, zlist_org, zlist_dyn, figpath, savename)
            %}
        ADC_breast_dyn = reshape(ADC_breast_dyn, [ysize, xsize, zsize]);
        
        
        disp('Align result of DWI mask ...')
        DWImask_breast_org_list = reshape(DWImask_breast(:,:,zidx_Dynamic2DWI_1:zidx_Dynamic2DWI_end),...
                                          [NumPix, zidx_Dynamic2DWI_end - zidx_Dynamic2DWI_1 + 1]);
        DWImask_breast_dyn = zeros(ysize*xsize, zsize);
        for p = 1:NumPix
            value_org = DWImask_breast_org_list(p,:)';
            value_dyn = interp1(Zlist_dwi_org, value_org, Zlist_dwi_dyn);
            DWImask_breast_dyn(p, zidx_DWI2Dynamic_1:zidx_DWI2Dynamic_end) = value_dyn';
        end
        DWImask_breast_dyn = double(DWImask_breast_dyn > 0.5);
            %%%%%%%% plot (optional)
            %{
            value_org = ADC_breast_org_list(CentralPix,:)';
            value_dyn = ADC_breast_dyn(CentralPix,max([zidx_DWI2Dynamic_1-5, 1]):min([zidx_DWI2Dynamic_end+5, zsize]))';
            zlist_org = Zlist_dwi_org;
            zlist_dyn = Zlist_dce(max([zidx_DWI2Dynamic_1-5, 1]):min([zidx_DWI2Dynamic_end+5, zsize]));
            savename = [CaseNum, '_Slicealignment_t1m-ADC_central'];
            Myplot_SliceResample_Compare(value_org, value_dyn, zlist_org, zlist_dyn, figpath, savename)
            %}
        DWImask_breast_dyn = reshape(DWImask_breast_dyn, [ysize, xsize, zsize]);
        
        
    
            %%%%% DWI image
        disp('Align result of DWI ...')
        dwi_breast_dyn = zeros(ysize, xsize, zsize, size(dwi_breast,4));
        for f = 1:size(dwi_breast,4)
            dwi_breast_org_list = reshape(dwi_breast(:,:,zidx_Dynamic2DWI_1:zidx_Dynamic2DWI_end,f),...
                                          [NumPix, zidx_Dynamic2DWI_end - zidx_Dynamic2DWI_1 + 1]);
            dwi_breast_dyn_f = zeros(ysize*xsize, zsize);
            for p = 1:NumPix
                value_org = dwi_breast_org_list(p,:)';
                value_dyn = interp1(Zlist_dwi_org, value_org, Zlist_dwi_dyn);
                dwi_breast_dyn_f(p, zidx_DWI2Dynamic_1:zidx_DWI2Dynamic_end) = value_dyn';
            end
                %%%%%%%% plot (optional)
                %{
                value_org = dwi_breast_org_list(CentralPix,:)';
                value_dyn = dwi_breast_dyn_f(CentralPix,max([zidx_DWI2Dynamic_1-5, 1]):min([zidx_DWI2Dynamic_end+5, zsize]))';
                zlist_org = Zlist_dwi_org;
                zlist_dyn = Zlist_dce(max([zidx_DWI2Dynamic_1-5, 1]):min([zidx_DWI2Dynamic_end+5, zsize]));
                savename = [CaseNum, '_Slicealignment_t1m-dwi_central_frame',num2str(f)];
                Myplot_SliceResample_Compare(value_org, value_dyn, zlist_org, zlist_dyn, figpath, savename)
                %}
            dwi_breast_dyn_f = reshape(dwi_breast_dyn_f, [ysize, xsize, zsize]);
            dwi_breast_dyn(:,:,:,f) = dwi_breast_dyn_f;
        end

        %close all  


        %% saving aligned datasets
        pos_x = Xlist_dce(xidx_DWI2Dynamic_1);
        pos_y = Ylist_dce(yidx_DWI2Dynamic_1);
        pos_z = Zlist_dce;

        OrigSetName_list  ={'DWI','DCE','ParaMap_ADC', 'ParaMap_DWImask_whole', 'ParaMap_DCEmask'};                         %% 'T1Mapping',  
        InputDataName_list={'dwi_breast_dyn','dce_breast','ADC_breast_dyn', 'DWImask_breast_dyn', 'DCEmask_breast'};    %% 't1m_breast',
        AlignSetName_list ={'DWI_align','DCE_align','ParaMap_ADC_align', 'ParaMap_DWImask_align', 'ParaMap_DCEmask_align'};       %% 'T1Mapping_align',

        for seriesID = 1:length(OrigSetName_list)
            OrigSetName   = OrigSetName_list{seriesID};
            InputDataName = InputDataName_list{seriesID};
            AlignSetName  = AlignSetName_list{seriesID};

            if exist(OrigSetName,'var')
                eval(['AlignSet  = ',OrigSetName,';'])
                eval(['InputData = ',InputDataName,';'])
                AlignSet.Data = InputData;
                AlignSet.pos  = repmat([pos_x*ones(size(pos_z)), pos_y*ones(size(pos_z)), pos_z],...
                                       [1,1,size(InputData,4)]);
                AlignSet.ori  = repmat(AlignSet.ori(1,:,1,1),...
                                       [zsize,1,size(InputData, 4)]);
                AlignSet.slloc= repmat(Zlist_dce, [1,1,size(InputData, 4)]);
                AlignSet.resolution = DCE.resolution;
                AlignSet.row        = ysize;
                AlignSet.columns    = xsize;
                AlignSet.slicethickness = DCE.slicethickness;

                eval([AlignSetName, ' = AlignSet;'])
                
                savename = [CaseNum, '_',AlignSetName,'.mat'];
                save([datapath, savename], AlignSetName,'-v7.3')
            end
        end

            [~, ~, Zlist_dce, EdgeXY_dce, ~, ~] = LocationExtract_ISPY2(DCE_align,flipped); 
            [~, ~, Zlist_dwi, EdgeXY_dwi, ~, ~] = LocationExtract_ISPY2(DWI_align,flipped);
        

            %%% plots (optional): FOV & Slice ovelapping
            %
            
            Myplot_FOVs_Overlapping  (EdgeXY_dce, EdgeXY_dwi, NaN, NaN, figpath, [CaseNum, 'align']) %%EdgeXY_t1m,
            Myplot_Slices_Overlapping(Zlist_dce,  Zlist_dwi,  NaN, NaN, figpath, [CaseNum, 'align']) %%Zlist_t1m,  
            %}
            
            %%%% plots (optional): central slices visualization & tranverse aligment comparison
            %
            DataSets = {dce_breast, dwi_breast_dyn, ADC_breast_dyn, DWImask_breast_dyn, DCEmask_breast};    %% t1m_breast, 
            Names    = {'DCE_align','DWI_align',    'ADC_align',    'DWImask_align',    'ParaMap_DCEmask'};                 %% 'T1Mapping', 
            
            [~, slc_dwi] = min(abs(Zlist_dwi - zloc));
            [~, slc_dce] = min(abs(Zlist_dce - zloc));
            Slices= [slc_dce, slc_dwi,  slc_dwi, slc_dwi, slc_dce];               %% slc_t1m,     

            Myplot_CentralSlices_Visual (CaseNum, figpath, DataSets, Names, Slices)
            Myplot_CentralSlices_Compare(CaseNum, figpath, DataSets, Names, Slices)
            %}

    end
beep;
end


%% visualization functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOV & Slice overlapping %%%%%%%%%%%%%%
function Myplot_FOVs_Overlapping(EdgeXY_dce, EdgeXY_dwi, EdgeXY_dti, EdgeXY_cest, figpath, CaseNum) %EdgeXY_t1m,

    figure
    hold on
    Names = {};
    plot(EdgeXY_dce(1:2,1), EdgeXY_dce(1:2,2),     'r-','LineWidth',4);
    Names{end+1} = 'DCE';
    
    plot(EdgeXY_dwi(1:2,1), EdgeXY_dwi(1:2,2),     'b-','LineWidth',4);
    Names{end+1} = 'DWI';
    
    plot(EdgeXY_dce(2:3,1), EdgeXY_dce(2:3,2),     'r-','LineWidth',4)
    plot(EdgeXY_dce(3:4,1), EdgeXY_dce(3:4,2),     'r-','LineWidth',4)
    plot(EdgeXY_dce([4,1],1), EdgeXY_dce([4,1],2), 'r-','LineWidth',4)
    
% %     plot(EdgeXY_t1m(1:2,1), EdgeXY_t1m(1:2,2),     'g-','LineWidth',1);
% %     plot(EdgeXY_t1m(2:3,1), EdgeXY_t1m(2:3,2),     'g-','LineWidth',1)
% %     plot(EdgeXY_t1m(3:4,1), EdgeXY_t1m(3:4,2),     'g-','LineWidth',1)
% %     plot(EdgeXY_t1m([4,1],1), EdgeXY_t1m([4,1],2), 'g-','LineWidth',1)
% %     Names{end+1} = 'T1M';

    
    plot(EdgeXY_dwi(2:3,1), EdgeXY_dwi(2:3,2),     'b-','LineWidth',4)
    plot(EdgeXY_dwi(3:4,1), EdgeXY_dwi(3:4,2),     'b-','LineWidth',4)
    plot(EdgeXY_dwi([4,1],1), EdgeXY_dwi([4,1],2), 'b-','LineWidth',4)
    

% %     if ~isnan(EdgeXY_cest(1))
% %         plot(EdgeXY_cest(1:2,1),EdgeXY_cest(1:2,2),    'm-','LineWidth',2);
% %         plot(EdgeXY_cest(2:3,1), EdgeXY_cest(2:3,2),   'm-','LineWidth',2)
% %         plot(EdgeXY_cest(3:4,1), EdgeXY_cest(3:4,2),   'm-','LineWidth',2)
% %         plot(EdgeXY_cest([4,1],1),EdgeXY_cest([4,1],2),'m-','LineWidth',2)
% %         Names{end+1} = 'CEST';
% %     end
% %         
% %     if ~isnan(EdgeXY_dti(1))
% %         plot(EdgeXY_dti(1:2,1), EdgeXY_dti(1:2,2),     'c-','LineWidth',1);
% %         plot(EdgeXY_dti(2:3,1), EdgeXY_dti(2:3,2),     'c-','LineWidth',1)
% %         plot(EdgeXY_dti(3:4,1), EdgeXY_dti(3:4,2),     'c-','LineWidth',1)
% %         plot(EdgeXY_dti([4,1],1), EdgeXY_dti([4,1],2), 'c-','LineWidth',1)
% %         Names{end+1} = 'DTI';
% %     end
    

    legend(Names,'Location','northeastoutside')    
        
    xlabel('x-coordinates [mm]')
    ylabel('y-coordinates [mm]')
%     xlim([-200 200])
%     ylim([-200 200])

    axis equal
    ax = gca;
    ax.FontSize = 22;
    ax.FontName = 'Times New Roman';
    set(ax, 'YDir','reverse')

    %% save figure
    if isempty(figpath)
        figpath = pwd;
    end
    
    savename = [CaseNum, '_FOVs'];
    saveas(gcf, [figpath, savename, '.jpg'])
    
end

function Myplot_Slices_Overlapping(Zlist_dce, Zlist_dwi, Zlist_dti, Zlist_cest, figpath, CaseNum) %Zlist_t1m,

    figure
    hold on
    Names = {};
    plot(1*ones(size(Zlist_dce)), Zlist_dce, 'r*-','LineWidth',1)
    Names{end+1} = 'DCE';
    
    plot(2*ones(size(Zlist_dwi)), Zlist_dwi, 'b*-','LineWidth',1)
    Names{end+1} = 'DWI';
    
    if ~isnan(Zlist_cest(1))
        plot(3*ones(size(Zlist_cest)),Zlist_cest,'m*-','LineWidth',1)
        Names{end+1} = 'CEST';
    end
    if ~isnan(Zlist_dti(1))
        plot(4*ones(size(Zlist_dti)), Zlist_dti, 'c*-','LineWidth',1)
        Names{end+1} = 'DTI';
    end
    

    ylabel('slices location [mm]')
    xlim([0.5 0.5+length(Names)])
    xticks(1:length(Names))
    xticklabels(Names)
    
    ax = gca;
    ax.FontSize = 22;
    ax.FontName = 'Times New Roman';

    
    %% save figure
    if isempty(figpath)
        figpath = pwd;
    end
    
%     savename = [CaseNum, '_Slices'];
%     saveas(gcf, [figpath, savename, '.fig'])
    savename = [CaseNum, '_Slices'];
    saveas(gcf, [figpath, savename, '.jpg'])
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% image visualization %%%%%%%%%%%%%%%%%
function Myplot_CentralSlices_Visual(CaseNum, figpath, DataSets, Names, Slices)
    for idx = 1:length(DataSets)
        Data = DataSets{idx};
        slc  = Slices(idx);
        vname= Names{idx};
        
        %
        figure
        for f = 1:size(Data, 4)
            img = Data(:,:,slc,f);
            high = max(img(:));
            low  = min(img(:));
            imshow(img,[low, high])
            savename = [CaseNum, '_',vname,'_central_frame', num2str(f)];
            saveas(gcf, [figpath, savename, '.jpg'])
        end
    end
    %close all
end

function Myplot_CentralSlices_Compare(CaseNum, figpath, DataSets, Names, Slices)
    figure
    Data = DataSets{1};
    img1 = Data(:,:,Slices(1),end);
    
    for idx = 2:length(DataSets)
        Data = DataSets{idx};
        slc  = Slices(idx);
        vname = Names{idx}
        
        f = size(Data, 4);
        if strcmp(vname, 'CESTcoef')
            f = 60;
        end
        img2 = Data(:,:,slc,f);
        imshowpair(img1, img2)
        savename = [CaseNum, '_FOValignment_',Names{1},'-',vname,'_central'];
        saveas(gcf, [figpath, savename, '.jpg'])
    end
    %close all

end
        
function Myplot_SliceResample_Compare(value_org, value_dyn, zlist_org, zlist_dyn, figpath, savename)
    figure
    plot(zlist_org, value_org, 'b*-', 'LineWidth',2, 'MarkerSize', 10)
    hold on
    plot(zlist_dyn, value_dyn, 'r*',  'LineWidth',2)
    legend('original slice', 'aligned slice')
    xlabel('z-coordinate of slices [mm]')
    xlim([zlist_dyn(end), zlist_dyn(1)])
    ylabel('signal [AU]')
    ax = gca;
    ax.FontSize = 22;
    ax.FontName = 'Times New Roman';

    saveas(gcf, [figpath, savename, '.jpg'])
    %close all
end

