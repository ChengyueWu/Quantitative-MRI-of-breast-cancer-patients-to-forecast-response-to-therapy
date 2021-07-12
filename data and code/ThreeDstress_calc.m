%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Function to calculate von mises stress 

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [stress_vm] = ThreeDstress_calc(U_disp,M_loc_m,bcf,dims,G,v)
% Inputs:
%      U_disp  :     displacement due to tumor cells double  (sy,sx,sz)    
%      M_loc_m :     Coordinates of Md matrix values double  <(sy,sx,sz)x3
%      bcf     :     boundary conditions             double  (sy,sx,sz)
%      dims    :     dx_TC, dy_TC, dz_TC             double  (1x3)
%      G       :     shear modulus                   double  (sy,sx,sz)
%      v       :     poisson's ratio                 double  (1x1)
% Outputs:
%      stress_vm:   calculated von mises stress      double  (sy,sx,sz) 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dx_TC = dims(1);
dy_TC = dims(2);
dz_TC = dims(3);

[sy sx sz] = size(G);
stress = zeros(sy,sx,sz,3);
U = zeros(sy,sx,sz,3);
count = size(M_loc_m,1);
    for j = 1:count
        y = M_loc_m(j,1);
        x = M_loc_m(j,2);
        z = M_loc_m(j,3);

        U(y,x,z,1) = U_disp(j);
        U(y,x,z,2) = U_disp(j+count);
        U(y,x,z,3) = U_disp(j+2*count);
    end
    U = -U;
    U = U+0;
   
    [strain_x,strain_xy,strain_xz] = Diffy3D(dx_TC,dy_TC,dz_TC,bcf,U(:,:,:,1));

    [strain_yx,strain_y,strain_yz] = Diffy3D(dx_TC,dy_TC,dz_TC,bcf,U(:,:,:,2));

    [strain_zx,strain_zy,strain_z] = Diffy3D(dx_TC,dy_TC,dz_TC,bcf,U(:,:,:,3));

    E = G.*(2*(1+v));
 
    lame = v*E;
    lame = lame/((1+v)*(1-2*v));
    
    lame = (2*G*v)/(1-2*v);
    
    stress(:,:,:,1) = lame.*(strain_x+strain_y+strain_z)+2*G.*strain_x;
    stress(:,:,:,2) = lame.*(strain_x+strain_y+strain_z)+2*G.*strain_y;
    stress(:,:,:,3) = lame.*(strain_x+strain_y+strain_z)+2*G.*strain_z;

    shear_xy = 0.5*G.*(abs(strain_xy) +abs(strain_yx));
    shear_xz = 0.5*G.*(abs(strain_xz) + abs(strain_zx));
    shear_yz = 0.5*G.*(abs(strain_yz) + abs(strain_zy));
    shear = 6*(shear_xy.^2+shear_xz.^2+shear_yz.^2);  %shear = 0;
   
    stress_vm = sqrt(0.5*((stress(:,:,:,1)-stress(:,:,:,2)).^2 + (stress(:,:,:,1)-stress(:,:,:,3)).^2 + (stress(:,:,:,2)-stress(:,:,:,3)).^2+shear));
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end of file