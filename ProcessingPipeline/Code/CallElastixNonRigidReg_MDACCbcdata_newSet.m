function CallElastixNonRigidReg_MDACCbcdata_newSet(subname,Parentpath,Elastixpath,target_list,Visit_num_f)
target_list = target_list - 1;
Visit_num_f = Visit_num_f - 1;
img_f_names = {[subname 'T',num2str(Visit_num_f),'_Init_hs.mhd'],... % image for initial reg
               [subname 'T',num2str(Visit_num_f),'_hs.mhd'],...      % intensity image
               };        
           
for Visit_num_m = target_list

    savename = [subname,'T',num2str(Visit_num_f),'T',num2str(Visit_num_m)]

    inpath  = [Parentpath, subname,'/'];
    outpath = [inpath subname 'T' num2str(Visit_num_m) 'regtoT' num2str(Visit_num_f)' '/'];
    if isfolder(outpath)==0
        mkdir(outpath)
    end
    
    img_m_names = {[subname 'T',num2str(Visit_num_m),'_Init_hs.mhd'],...
                   [subname 'T',num2str(Visit_num_m),'_hs.mhd'],...
                   };

% %     parafile_names = {'Rigid.txt',['BSpline_Rigid_PenaltyV' num2str(Visit_num_m) '.txt']};
    parafile_name = ['BSpline_Rigid_PenaltyT', num2str(Visit_num_m), '_',subname '.txt'];

    % 1. first, register the image with rigid registration as initial alignment
    idx = 1;
    img_m     = img_m_names{1};
    img_f     = img_f_names{1};
    
    info       = mhd_read_header([inpath, img_m]);
    InitMask_m = mhd_read_volume(info);
    info       = mhd_read_header([inpath, img_f]);
    InitMask_f = mhd_read_volume(info);
    
    [optimizer, metric] = imregconfig('multimodal');
    metric.NumberOfSpatialSamples = 5000;
    optimizer.InitialRadius = 0.0003;
    optimizer.GrowthFactor  = 1.0001;
    optimizer.MaximumIterations = 500;
    
    [InitMask_reg, transform3D, Ref] = RigidRegister_3D(InitMask_m, InitMask_f, [], [1,1,1]);
    DeformFields{idx}.D   = transform3D;
    DeformFields{idx}.Ref = Ref;
    
    img_apply = img_m_names{2};
    info = mhd_read_header([inpath, img_apply]);
    Target_org = mhd_read_volume(info);
    Target_reg = imwarp(Target_org, Ref, transform3D, 'OutputView',Ref);

    Convert_MAT2MHD_ForElastix(Target_reg, outpath, 'outputname', 'result')
    
    
% %     idx = 1;
% %     img_m     = img_m_names{2};
% %     img_apply = img_m_names{2};
% %     command = ['cd ' inpath,';',... go to the folder containing images and parameter files
% %                'export PATH=' Elastixpath 'bin/:$PATH;', ...
% %                'export DYLD_LIBRARY_PATH=' Elastixpath 'lib:$DYLD_LIBRARY_PATH;',...
% %                Elastixpath 'bin/elastix',...
% %                ' -f ', img_f_names{2},...
% %                ' -m ', img_m,...
% %                ' -out ',outpath,...
% %                ' -p ',Parentpath, parafile_names{1},';']; 
% %     [flag1,flag2] = system(command,'-echo')
% % 
% %     % 1.1 generate the deformation field of the initial alignment & apply to the moving image
% %     %%% !!! note that before generating the deformation field from the binary
% %     %%% masks alignment and applying to intensity images, interporator order in 
% %     %%% the transform parameter file (obtained from step 1) need to be changed
% %     disp('Change B-spline interporator order: from 0 to 3 ...')
% %     pause(1)
% %     fid = fopen([outpath 'TransformParameters.0.txt']);
% %     C = textscan(fid,'%s','delimiter','\n');
% %     C = C{1};
% %     fclose(fid);
% %     fid = fopen([outpath 'TransformParameters.0.txt'], 'w');
% %     for k = 1:numel(C)
% %         command = C{k};
% %         if regexp(command,'FinalBSplineInterpolationOrder')
% %             command = '(FinalBSplineInterpolationOrder 3)';
% %         end
% %         fprintf(fid, '%s\r\n', command);
% %     end
% %     fclose(fid);
% % 
% %     command = ['cd ' inpath,';',... go to the folder containing images and parameter files
% %                'export PATH=' Elastixpath 'bin/:$PATH;', ...
% %                'export DYLD_LIBRARY_PATH=' Elastixpath 'lib:$DYLD_LIBRARY_PATH;',...
% %                Elastixpath 'bin/transformix -in ',img_apply,' -def all -out ',outpath,' -tp ',outpath,'TransformParameters.0.txt']; %% calling elastix
% %     [~,~] = system(command,'-echo');
% %     [DeformFields{idx},~] = read_mhd([outpath, 'deformationField.mhd']);

    % 2. second, register with B-Spline + rigidity penalty, and generate the second deformation field
    idx = 2;
    img_m = 'result.mhd'; %% 'result.mhd': moving image registered by the initial rigid alignment
    img_f = img_f_names{2};
    
    command = ['cd '    outpath,';',... go to the folder containing images and parameter files
               'export PATH=' Elastixpath 'bin/:$PATH;', ...
               'export DYLD_LIBRARY_PATH=' Elastixpath 'lib:$DYLD_LIBRARY_PATH;',...
               Elastixpath 'bin/elastix',...
               ' -f ',  inpath,  img_f,...
               ' -m ',  outpath, img_m,... img_m_names{2}
               ' -out ',outpath,... 
               ' -p ',  inpath, parafile_name,';',... parafile_names{2}
               Elastixpath 'bin/transformix -jac all -def all -out ',outpath,' -tp ',outpath,'TransformParameters.0.txt']; %% calling elastix

    [~,~] = system(command,'-echo');

    [DeformFields{idx},~] = read_mhd([outpath, 'deformationField.mhd']);
    [img_final, ~]        = read_mhd([outpath, 'result.0.mhd']);
    [spatialJacobian, ~]  = read_mhd([outpath, 'spatialJacobian.mhd']);

    save([outpath savename '.mat'],'img_final','DeformFields','spatialJacobian', '-v7.3');
    
end
