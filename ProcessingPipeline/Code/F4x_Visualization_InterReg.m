clc
clear 
close all

Scriptpath = '/Users/chengyuewu/Desktop/study_work/Lab/MDACC_ClinicalPipeline/MDACC_patient_processing/Codes/Clinical-Breast-Cancer-Model-Development-master/MDACC_DataProcessing_Reports/Scripts/';
matpath  = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';

Patient_list = {
'104268'...
};

dx = 0.7031; dy = 0.7031; dz = 2;
          
%%
for idx = 1:length(Patient_list)
    
    subname = Patient_list{idx};
    subname
    
    % cd(Scriptpath)
    inpath  = [matpath, subname, '/'];
    figpath = [matpath, subname, '/figures/'];
    if ~isfolder(figpath)
        mkdir(figpath)
    end
    
    T0_org = load([inpath, subname, 'T0_MovingToBeRegTo_T1.mat']);
    T0_reg = load([inpath, subname, 'T0_Reg.mat']);
    T2_org = load([inpath, subname, 'T2_MovingToBeRegTo_T1.mat']);
    T2_reg = load([inpath, subname, 'T2_Reg.mat']);
    T1     = load([inpath, subname, 'T1_Reg.mat']);
    
    ctr_s = round(median(find(sum(sum(T1.roi)))));
    close all

    %% check 3D alignment of tumor masks
    %
    h1 = figure;
    
    testvol1 = logical(T0_org.roi);
    col=[0 0 0.8];
    hiso1 = patch(isosurface(testvol1,0),'FaceColor',col,'EdgeColor','none');
    alpha(0.3);
    hold on
    testvol1 = logical(T1.roi);
    col=[0 0.8 0];
    hiso1 = patch(isosurface(testvol1,0),'FaceColor',col,'EdgeColor','none');
    alpha(0.3);
    hold on
    testvol1 = logical(T2_org.roi);
    col=[0.8 0 0];
    hiso1 = patch(isosurface(testvol1,0),'FaceColor',col,'EdgeColor','none');
    alpha(0.3);
    
    lighting phong;
    camlight;
    view([5 -2 1])
    grid on
    axis equal tight
    legend({'Visit 1', 'Visit 2', 'Visit 3'},'location','northeast')
    set(gca,'DataAspectRatio',[1/dx 1/dy 1/dz])
    % 
    % saveas(h, [figpath subname, '_3Dtumors.fig'])
    % saveas(h, [figpath subname, '_3Dtumors.jpg'])

    
    %%%%%%%%
    h2 = figure;
    
    testvol1 = logical(T0_reg.roi);
    col=[0 0 0.8];
    hiso1 = patch(isosurface(testvol1,0),'FaceColor',col,'EdgeColor','none');
    alpha(0.3);
    hold on
    testvol1 = logical(T1.roi);
    col=[0 0.8 0];
    hiso1 = patch(isosurface(testvol1,0),'FaceColor',col,'EdgeColor','none');
    alpha(0.3);
    hold on
    testvol1 = logical(T2_reg.roi);
    col=[0.8 0 0];
    hiso1 = patch(isosurface(testvol1,0),'FaceColor',col,'EdgeColor','none');
    alpha(0.3);
    
    lighting phong;
    camlight;
    view([5 -2 1])
    grid on
    axis equal tight
    legend({'Visit 1', 'Visit 2', 'Visit 3'},'location','northeast')
    set(gca,'DataAspectRatio',[1/dx 1/dy 1/dz])
    % 
    % saveas(h, [figpath subname, '_3Dtumors.fig'])
    % saveas(h, [figpath subname, '_3Dtumors.jpg'])
    %}
    
    %% check central slices with tumor masks
    %
    %%%%%%%
    h1 = figure;
    set(h1, 'Unit','characters','Position', [10, 30, 230, 30]);
    
    subplot(1,3,1)
    imagesc(T0_org.avgdce(:,:,ctr_s))
    axis off equal tight
    tumoredge = edge(T0_org.roi(:,:,ctr_s));
    [yt, xt]       = find(tumoredge);
    hold on; plot(xt,yt,'r.')
    
    subplot(1,3,2)
    imagesc(T1.avgdce(:,:,ctr_s))
    axis off equal tight
    tumoredge = edge(T1.roi(:,:,ctr_s));
    [yt, xt]       = find(tumoredge);
    hold on; plot(xt,yt,'r.')
    
    subplot(1,3,3)
    imagesc(T2_org.avgdce(:,:,ctr_s))
    axis off equal tight
    tumoredge = edge(T2_org.roi(:,:,ctr_s));
    [yt, xt]       = find(tumoredge);
    hold on; plot(xt,yt,'r.')
    
    % saveas(h, [figpath subname, '_ctr_DCEtumors.jpg'])

    %%%%%%%
    h2 = figure;
    set(h2, 'Unit','characters','Position', [10, 30, 230, 30]);
    
    subplot(1,3,1)
    imagesc(T0_reg.avgdce(:,:,ctr_s))
    axis off equal tight
    tumoredge = edge(T0_reg.roi(:,:,ctr_s));
    [yt, xt]       = find(tumoredge);
    hold on; plot(xt,yt,'r.')
    
    subplot(1,3,2)
    imagesc(T1.avgdce(:,:,ctr_s))
    axis off equal tight
    tumoredge = edge(T1.roi(:,:,ctr_s));
    [yt, xt]       = find(tumoredge);
    hold on; plot(xt,yt,'r.')
    
    subplot(1,3,3)
    imagesc(T2_reg.avgdce(:,:,ctr_s))
    axis off equal tight
    tumoredge = edge(T2_reg.roi(:,:,ctr_s));
    [yt, xt]       = find(tumoredge);
    hold on; plot(xt,yt,'r.')
    
    % saveas(h, [figpath subname, '_ctr_DCEtumors.jpg'])
    
    %}
    
    %% check deformation grid 
    for v = [0,2]
        InputFilename  = [subname, 'T',num2str(v),'_MovingToBeRegTo_T1.mat'];
        load([inpath, InputFilename])
        
        DeformFilename = [subname, 'T1T',num2str(v), '.mat'];
        load([inpath, subname, 'T',num2str(v),'regtoT1/', DeformFilename])
        
        ctr = round(median(find(sum(sum(roi_orig)))));
        TumorEdge        = edge(roi_orig(:,:,ctr));
        [yTE, xTE]       = find(TumorEdge);
        
        eval(['roi_df = T',num2str(v),'_reg.roi;'])
        ctr_df = round(median(find(sum(sum(roi_df)))));
        TumorEdge        = edge(roi_df(:,:,ctr_df));
        [yTE_df, xTE_df]       = find(TumorEdge);
        
        test_grid = zeros(img_final.size');
            for z = 1:img_final.size(3)
                test_grid(1:5:img_final.size(1), :, z) = z; %img_final.size(3)-z*2;
                test_grid(:,1:5:img_final.size(2), z)  = z; %img_final.size(3)-z*2;
            end

            test_grid_df = test_grid;

                %%%%%% initial alignment
                D   = DeformFields{1}.D;
                Ref = DeformFields{1}.Ref;
                test_grid_df = imwarp(test_grid_df,Ref, D, 'OutputView',Ref);

                %%%%%% non-rigid reg
                for df_idx = 2:length(DeformFields)
                    D = zeros(size(DeformFields{df_idx}.data,1),size(DeformFields{df_idx}.data,2),size(DeformFields{df_idx}.data,3),3);
                    D(:,:,:,1) = DeformFields{df_idx}.datay;
                    D(:,:,:,2) = DeformFields{df_idx}.datax;
                    D(:,:,:,3) = DeformFields{df_idx}.dataz;

                    test_grid_df = imwarp(test_grid_df, D);
                end

            %
            figure;
            subplot(1,2,1)
            imagesc(test_grid(:,:,ctr))
            colormap(gray)  %     colorbar
            title('original central tumor slice')
            axis off
            axis equal tight
            hold on
            plot(xTE,yTE,'.r')
            ax = gca;
            ax.FontSize = 20;
            ax.FontName = 'Times New Roman';

            subplot(1,2,2)
            imagesc(test_grid_df(:,:,ctr_s))
            colormap(gray)      %     colorbar
            title('deformed central tumor slice ')
            axis off
            axis equal tight
            hold on
            plot(xTE_df,yTE_df,'.r')
            ax = gca;
            ax.FontSize = 20;
            ax.FontName = 'Times New Roman';

            set(gcf,'Units','characters','Position',[1 10 100 30]);
        %     saveas(h1,[test_path subname '_' method, '_v',num2str(visit_num_f),'v',num2str(visit_num_m),'_testdf_ctrTUMOR.fig'])
        %     saveas(gcf,[figpath subname '_' method, '_v',num2str(visit_num_f),'v',num2str(v),'_testdf_ctrTUMOR.jpg'])
    end
    
end

           
