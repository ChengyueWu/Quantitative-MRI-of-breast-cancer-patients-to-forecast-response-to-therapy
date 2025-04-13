%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Function to calculate the mechanical coupling to the diffusion and call
%the forward finite difference method to calculate the number of tumor
%cells for each time step

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [TCs] = NtcFDmech3DUTLogisticOnly(TCs_in,dims,Ds,dt_TC,k,carcap,bcf,steps,G,DiffyG,v)
% Inputs:
%      TCs_in  :     Tumor Cells In                  double  (sy,sx,sz)
%      dims    :     dx_TC, dy_TC, dz_TC             double  (1x3)
%      dx_TC   :     grid spacing in x               double  (1x1)
%      dy_TC   :     grid spacing in y               double  (1x1)
%      dz_TC   :     grid spacing in z               double  (1x1)
%      Ds      :     Tumor Cell Diffusion            double  (sy,sx,sz)
%      dt_TC   :     time step                       double  (1x1)
%      k       :     growth rate                     double  (sy,sx,sz)
%      carcap  :     Carrying Capacity               double  (sy,sx,sz)
%      bcf     :     boundary conditions             double  (sy,sx,sz)
%      steps   :     time steps                      double  (1x1)
%      G       :     shear modulus                   double  (sy,sx,sz)
%      DiffyG  :     derivative of G                 double  (sy,sx,sz)
%      v       :     poisson's ratio                 double  (1x1)
% Outputs:
%      TCs     :     Tumor Cells out                 double(sy,sx,sz)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Defining the dimensions
dx_TC = dims(1); 
dy_TC = dims(2); 
dz_TC = dims(3);

[sy_TC, sx_TC, sz_TC]  = size(TCs_in);

%Holding matrices
TCs = zeros(sy_TC,sx_TC,sz_TC);
DiffyN = zeros(size(DiffyG));

%Defining fixed parameters for the linear isotropic equilibrium of the
%mechanical forces
lam_1_e = 2.5e-3;
lam_2_e = 2.5e-3;

%Only calculate the mechanical coupling changes every 20 time steps
mechstep = 20;
for inncount = [1,mechstep:mechstep:steps]
    if(inncount==1||mod(inncount,mechstep)==0)
        %Calculate the derivative
        [Diffyx,Diffyy,Diffyz] = Diffy3D(dx_TC,dy_TC,dz_TC,bcf,TCs_in);
        DiffyN(:,:,:,1) = Diffyx;
        DiffyN(:,:,:,2) = Diffyy;
        DiffyN(:,:,:,3) = Diffyz;
        %Call functions for deriving the mechanical stresses
        [F_list,Md,M_loc_m] = ...
                  ThreeDmeq_matrix_builder_opt(bcf,G,v,dims,DiffyG,DiffyN);
        setup.type = 'crout'; 
        setup.droptol = 1e-4; 
        setup.milu = 'row';        
        [svm] = ThreeDmech_opt_solver_bcg(Md,DiffyN,F_list,...
                                        (10^lam_1_e),M_loc_m,G,v,dims,bcf);            
        Dtc = Ds.*exp(-lam_2_e*svm); 
    end
    
    %Updating steps 
    if inncount == 1
        nextstep = mechstep-2;
    elseif (inncount + mechstep-1) < steps
        nextstep = mechstep-1;
    else
        nextstep = steps-inncount;
    end
    
    %Call the forward solver for the time step
    [TCs] = forwardsolveLogisticOnly(TCs_in,dt_TC,k,carcap,bcf,sy_TC,...
                     sx_TC,sz_TC,Dtc,dx_TC,dy_TC,dz_TC,inncount,nextstep);
    %Update tumor cell number             
    TCs_in = TCs;  
end %end of time step loop   
 
end %end of function

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end of file