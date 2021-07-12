%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Function to calculate von mises stress

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [stress_vm] = ThreeDmech_opt_solver_bcg(Md,Fgrad,F_list,mech_lam,M_loc_m,G,v,dims,bcf)
% Inputs:
%      Md      :     FD coeff matrix for mech        double  <(sy,sx,sz)
%      Fgrad   :     derivative of tumor cells       double  (sy,sx,sz)
%      F_list  :     Gradient of tumor cells         double  (1x<sy*sx*sz*3)
%      mech_lam:     coupling constant               double  (1x1)
%      M_loc_m :     Coordinates of Md matrix values double  <(sy,sx,sz)x3
%      G       :     shear modulus                   double  (sy,sx,sz)
%      v       :     poisson's ratio                 double  (1x1)
%      dims    :     dx_TC, dy_TC, dz_TC             double  (1x3)
%      bcf     :     boundary conditions             double  (sy,sx,sz)
% Outputs:
%      stress_vm:   calculated von mises stress      double  (sy,sx,sz) 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Splitting up the derivatives
Fx = Fgrad(:,:,:,1); 
Fy = Fgrad(:,:,:,2); 
Fz = Fgrad(:,:,:,3); 

%
Fv = Fx(F_list);
Fv = [Fv;Fy(F_list)];
Fv = [Fv;Fz(F_list)];

%Calculate the displacement from tumor cells (allowing extra outputs
%eliminates output in the command window)
[U_disp flag relres iters resvec]= bicgstab(Md,-1*mech_lam*Fv,10^-2.2,1e4);

%Calculate the resulting von mises stress
[stress_vm] = ThreeDstress_calc(U_disp,M_loc_m,bcf,dims,G,v);
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end of file