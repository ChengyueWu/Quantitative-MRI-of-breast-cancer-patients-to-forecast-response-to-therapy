function CallElastixNonRigidReg_ISPY2(subname,Parentpath,Elastixpath,target_list,Visit_num_f)
target_list = target_list - 1;
Visit_num_f = Visit_num_f - 1;


img_f_names = {[subname 'T',num2str(Visit_num_f),'_Init_hs.mhd'],... % image for initial reg
               [subname 'T',num2str(Visit_num_f),'_hs.mhd'],...      % intensity image
               };        
           
for Visit_num_m = target_list

    savename = [subname,'T',num2str(Visit_num_f),'T',num2str(Visit_num_m)]

    inpath = [Parentpath, subname,'/'];
    % inpath = ['/Volumes/SanDiskSSD/ISPY2/Data/Processed/', subname,'/IntervisitRegistered/RegFiles/'];

    outpath = [inpath subname 'T' num2str(Visit_num_m) 'regtoT' num2str(Visit_num_f)' '/'];
    if isfolder(outpath)==0
        mkdir(outpath)
    end
    
    img_m_names = {[subname 'T',num2str(Visit_num_m),'_Init_hs.mhd'],...
                   [subname 'T',num2str(Visit_num_m),'_hs.mhd'],...
                   };

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
    
    
    mask_apply = [subname 'T',num2str(Visit_num_m),'_TumorMask_hs'];
    info = mhd_read_header([inpath, mask_apply, '.mhd']);
    mask_orig  = mhd_read_volume(info);
    mask_reg   = imwarp(mask_orig, Ref, transform3D, 'OutputView',Ref);
    mask_reg   = double(mask_reg > 0.5);
    Convert_MAT2MHD_ForElastix(mask_reg, inpath, 'outputname', mask_apply)
    
    
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
