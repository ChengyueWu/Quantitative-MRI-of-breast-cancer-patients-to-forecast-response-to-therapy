%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Function to calculate diffusion matrix

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [F_list,Md,M_loc_m] = ThreeDmeq_matrix_builder_opt(bcf,G,v,dims,Ggrad,Fgrad)
% Inputs:
%      bcf     :     boundary conditions            double  (sy,sx,sz)
%      G       :     shear modulus                  double  (sy,sx,sz)
%      v       :     poisson's ratio                double  (1x1)
%      dims    :     dx_TC, dy_TC, dz_TC            double  (1x3)
%      Ggrad   :     derivative of G                double  (sy,sx,sz)
%      Fgrad   :     derivative of tumor cells      double  (sy,sx,sz)
% Outputs:
%      F_list  :     Gradient of tumor cells        double  (1x<sy*sx*sz*3)
%      Md      :     FD coeff matrix for mech       double  <(sy,sx,sz)
%      M_loc_m :     Coordinates of Md matrix values double  <(sy,sx,sz)x3
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fym = [1 11 21 101 111 121 201 211 221 301 311 321];
fyp = [2 12 22 102 112 122 202 212 222 302 312 322];
fy  = [fym fyp];

fxm = [20 21 22 120 121 122 220 221 222 320 321 322];
fxp = [10 11 12 110 111 112 210 211 212 310 311 312];
fx  = [fxm fxp];

fzm = [200 201 202 205 210 211 212 220 221 222 320 321 322];
fzp = [100 101 102 105 110 111 112 120 121 122 300 301 302 305];
fz  = [fzm fzp];
                
                
Gx = Ggrad(:,:,:,1);
Gy = Ggrad(:,:,:,2);
Gz = Ggrad(:,:,:,3);
Fx = Fgrad(:,:,:,1);
Fy = Fgrad(:,:,:,2);
Fz = Fgrad(:,:,:,3);

[sy sx sz] = size(bcf);

bct = ones(size(bcf));

Mx = zeros(sum(bcf(:)==5)+sum(bct(:)==2),19);
My = Mx; Mz = Mx;
Fx_v = zeros(size(Mx,1),1);
F_list = zeros(size(Mx,1),1);
Fy_v = Fx_v;
Fz_v = Fy_v;

Ux_bc = zeros(size(Mx,1),1);
Uy_bc = Ux_bc;
Uz_bc = Uy_bc;

dx = dims(1); dy = dims(2); dz = dims(3);

count = 0;
M_loc = zeros(sy,sx,sz);
k = (1/(1*2*v));
j = 0;
xtype = [10 11 12 20 21 22 110 111 112 120 121 122 210 211 212 220 221 222];
ytype = [1  2  11 12 21 22 101 102 111 112 121 122 201 202 211 212 221 222];
ztype = [100 101 102 105 110 111 112 120 121 122 200 201 202 205 210 211 212 220 221 222 300 301 302 305 310 311 312 320 321 322];
for z = 1:sz
    for x = 1:sx 
        for y = 1:sy  
                    
            f = bcf(y,x,z);
            ft = bct(y,x,z);
            j = j+1;
            
            if f == 5
                count = count +1;
                M_loc(j) = count;
                Ux_bc(count,1) = 0;
                Uy_bc(count,1) = 0;
                Uz_bc(count,1) = 0;
                
                % ---- Mx Build ------
                Mx(count,1) = -2*G(j)*((1/dx^2)+(1/dy^2)+(1/dz^2)+(k/dx^2));
                
                Mx(count,2) = (Gx(j)/(2*dx))*(1+k)+(G(j)/dx^2)*(1+k);
                nb = bcf(y,x+1,z);
                if nb ==0
                    bctype  = bct(y,x+1,z);
                    if bctype == 1
                        Ux_bc(count,1) = Ux_bc(count,1);
                        Mx(count,2) = 0;
                    end
                elseif sum(nb == xtype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,2) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end
                
                               
                Mx(count,3) = -(Gx(j)/(2*dx))*(1+k)+(G(j)/dx^2)*(1+k);
                nb = bcf(y,x-1,z);
                if nb ==0
                    bctype  = bct(y,x-1,z);
                    if bctype == 1
                        Ux_bc(count,1) = Ux_bc(count,1);
                        Mx(count,3) = 0;
                    end
                elseif sum(nb == xtype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,3) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end
                
                
                Mx(count,4) = (Gy(j)/(2*dy)) + G(j)/dy^2;
                nb = bcf(y+1,x,z);
                if nb ==0
                    bctype  = bct(y+1,x,z);
                    if bctype == 1
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,4) = 0;
                    end
                elseif sum(nb == xtype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,4) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end              
                               
                Mx(count,5) = (-Gy(j)/(2*dy)) + G(j)/dy^2;
                nb = bcf(y-1,x,z);
                if nb ==0
                   bctype = bct(y-1,x,z);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,5) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,5) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end                
                               
                Mx(count,6) = (Gz(j)/(2*dz)) + G(j)/dz^2;
                nb = bcf(y,x,z+1);
                if nb ==0
                   bctype = bct(y,x,z+1);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,6) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,6) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end          
                
                Mx(count,7) = (-Gz(j)/(2*dz)) + G(j)/dz^2;
                nb = bcf(y,x,z-1);
                if nb ==0
                   bctype = bct(y,x,z-1);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,7) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,7) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end          
                  
                Mx(count,8) = k*Gx(j)/(2*dy);
                nb = bcf(y+1,x,z);
                if nb ==0
                   bctype = bct(y+1,x,z);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,8) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,8) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end                
                
                Mx(count,9) = -k*Gx(j)/(2*dy);
                nb = bcf(y-1,x,z);
                if nb ==0
                   bctype = bct(y-1,x,z);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,9) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,9) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end       
                
                Mx(count,10) = k*G(j)/(4*dx*dy);
                nb = bcf(y+1,x+1,z);
                if nb ==0
                   bctype = bct(y+1,x+1,z);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,10) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,10) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end                
                
                Mx(count,11) = -k*G(j)/(4*dx*dy);
                nb = bcf(y-1,x+1,z);
                if nb ==0
                   bctype = bct(y-1,x+1,z);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,11) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,11) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end       
                
                Mx(count,12) = -k*G(j)/(4*dx*dy);
                nb = bcf(y+1,x-1,z);
                if nb ==0
                   bctype = bct(y+1,x-1,z);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,12) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,12) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end                
                
                Mx(count,13) = k*G(j)/(4*dx*dy);
                nb = bcf(y-1,x-1,z);
                if nb ==0
                   bctype = bct(y-1,x,z);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,13) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,13) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end       
                
                Mx(count,14) = Gx(j)*k/(2*dz);
                nb = bcf(y,x,z+1);
                if nb ==0
                   bctype = bct(y,x,z+1);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,14) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,14) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end                
                
                Mx(count,15) = -Gx(j)*k/(2*dz);
                nb = bcf(y,x,z-1);
                if nb ==0
                   bctype = bct(y,x,z-1);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,15) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,15) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end       
                
                Mx(count,16) = k*G(j)/(4*dx*dz);
                nb = bcf(y,x+1,z+1);
                if nb ==0
                   bctype = bct(y,x+1,z+1);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,16) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,16) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end                
                
                Mx(count,17) = -k*G(j)/(4*dx*dz);
                nb = bcf(y,x+1,z-1);
                if nb ==0
                   bctype = bct(y,x+1,z-1);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,17) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,17) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end       
                
                Mx(count,18) = -k*G(j)/(4*dx*dz);
                nb = bcf(y,x-1,z+1);
                if nb ==0
                   bctype = bct(y,x-1,z+1);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,18) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,18) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end                
                
                Mx(count,19) = k*G(j)/(4*dx*dz);
                nb = bcf(y,x-1,z-1);
                if nb ==0
                   bctype = bct(y-1,x,z);
                   if bctype == 1 %dirchlett
                      Ux_bc(count,1) = Ux_bc(count,1);
                      Mx(count,19) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Ux_bc(count,1) = Ux_bc(count,1)+0;
                    Mx(count,19) = 0;
                else
                    Ux_bc(count,1) = Ux_bc(count,1) + 0;
                end
                
                Fx_v(count,1) = Fx(j);

                
  
                % ---- My Build ------
                My(count,1) = -2*G(j)*((1/dx^2)+(1/dy^2)+(1/dz^2)+(k/dy^2));
                
                My(count,2) = (Gy(j)/(2*dy))*(1+k)+(G(j)/dy^2)*(1+k);
                nb = bcf(y+1,x,z);
                if nb ==0
                   bctype = bct(y+1,x,z);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,2) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,2) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end
                
                My(count,3) = -(Gy(j)/(2*dy))*(1+k)+(G(j)/dy^2)*(1+k);
                nb = bcf(y-1,x,z);
                if nb ==0
                   bctype = bct(y-1,x,z);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,3) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,3) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end     
                
                My(count,4) = (Gx(j)/(2*dx))+G(j)/dx^2;
                nb = bcf(y,x+1,z);
                if nb ==0
                   bctype = bct(y,x+1,z);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,4) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,4) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end                     
                
                My(count,5) = -(Gx(j)/(2*dx))+G(j)/dx^2;
                nb = bcf(y,x-1,z);
                if nb ==0
                   bctype = bct(y,x-1,z);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,5) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,5) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end     
                 
                My(count,6) = (Gz(j)/(2*dz))+G(j)/dz^2;
                nb = bcf(y,x,z+1);
                if nb ==0
                   bctype = bct(y,x,z+1);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,6) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,6) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end                     
                
                My(count,7) = -(Gz(j)/(2*dz))+G(j)/dz^2;
                nb = bcf(y,x,z-1);
                if nb ==0
                   bctype = bct(y,x,z-1);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,7) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,7) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end   
                
                My(count,8) = k*Gy(j)/(2*dx);
                nb = bcf(y,x+1,z);
                if nb ==0
                   bctype = bct(y,x+1,z);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,8) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,8) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end                     
                
                My(count,9) = -k*Gy(j)/(2*dx);
                nb = bcf(y,x-1,z);
                if nb ==0
                   bctype = bct(y,x-1,z);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,9) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,9) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end    
                
                My(count,10) = k*G(j)/(4*dx*dy);
                nb = bcf(y+1,x+1,z);
                if nb ==0
                   bctype = bct(y+1,x+1,z);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,10) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,10) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end                     
                
                My(count,11) = -k*G(j)/(4*dx*dy);
                nb = bcf(y+1,x-1,z);
                if nb ==0
                   bctype = bct(y+1,x-1,z);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,11) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,11) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end  
                
                My(count,12) = -k*G(j)/(4*dx*dy);
                nb = bcf(y-1,x+1,z);
                if nb ==0
                   bctype = bct(y-1,x+1,z);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,12) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,12) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end                     
                
                My(count,13) = k*G(j)/(4*dx*dy);
                nb = bcf(y-1,x-1,z);
                if nb ==0
                   bctype = bct(y-1,x-1,z);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,13) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,13) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end  
                
                My(count,14) = k*Gy(j)/(2*dz);
                nb = bcf(y,x,z+1);
                if nb ==0
                   bctype = bct(y,x,z+1);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,14) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,14) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end                     
                
                My(count,15) = -k*Gy(j)/(2*dz);
                nb = bcf(y,x,z-1);
                if nb ==0
                   bctype = bct(y,x,z-1);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,15) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,15) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end   
                
                My(count,16) = k*G(j)/(4*dz*dy);
                nb = bcf(y+1,x,z+1);
                if nb ==0
                   bctype = bct(y+1,x,z+1);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,16) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,16) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end                     
                
                My(count,17) = -k*G(j)/(4*dz*dy);
                nb = bcf(y+1,x,z-1);
                if nb ==0
                   bctype = bct(y+1,x,z-1);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,17) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,17) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end   
                
                My(count,18) = -k*G(j)/(4*dz*dy);
                nb = bcf(y-1,x,z+1);
                if nb ==0
                   bctype = bct(y-1,x,z+1);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,18) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,18) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end                     
                
                My(count,19) = k*G(j)/(4*dz*dy);
                nb = bcf(y-1,x,z-1);
                if nb ==0
                   bctype = bct(y-1,x,z-1);
                   if bctype == 1 %dirchlett
                      Uy_bc(count,1) = Uy_bc(count,1);
                      My(count,19) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uy_bc(count,1) = Uy_bc(count,1)+0;
                    My(count,19) = 0;
                else
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                end   
                
                Fy_v(count,1) = Fy(j);
                
  
                % ---- Mz Build ------
                Mz(count,1) = -2*G(j)*((1/dx^2)+(1/dy^2)+(1/dz^2)+(k/dz^2));
                
                Mz(count,2) = (G(j)/(dx^2))+Gx(j)/(2*dx);
                nb = bcf(y,x+1,z);
                if nb ==0
                   bctype = bct(y,x+1,z);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,2) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,2) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,3) = (G(j)/(dx^2))-Gx(j)/(2*dx);
                nb = bcf(y,x-1,z);
                if nb ==0
                   bctype = bct(y,x-1,z);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,3) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,3) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,4) = (G(j)/(dy^2))+Gy(j)/(2*dy);
                nb = bcf(y+1,x,z);
                if nb ==0
                   bctype = bct(y+1,x,z);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,4) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,4) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end                
                
                Mz(count,5) = (G(j)/(dy^2))-Gy(j)/(2*dy);
                nb = bcf(y-1,x,z);
                if nb ==0
                   bctype = bct(y-1,x,z);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,5) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,5) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end     
                
                Mz(count,6) = (Gz(j)/(2*dz))*(1+k)+(G(j)/dz^2)*(1+k);
                nb = bcf(y,x,z+1);
                if nb ==0
                   bctype = bct(y,x,z+1);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,6) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,6) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end   

                Mz(count,7) = -(Gz(j)/(2*dz))*(1+k)+(G(j)/dz^2)*(1+k);
                nb = bcf(y,x,z-1);
                if nb ==0
                   bctype = bct(y,x,z-1);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,7) = 0;
                   end
                elseif sum(nb == ztype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,7) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end   

                Mz(count,8) = k*Gz(j)/(2*dx);
                nb = bcf(y,x+1,z);
                if nb ==0
                   bctype = bct(y,x+1,z);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,8) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,8) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,9) = -k*Gz(j)/(2*dx);
                nb = bcf(y,x-1,z);
                if nb ==0
                   bctype = bct(y,x-1,z);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,9) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,9) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,10) = k*G(j)/(4*dz*dx);
                nb = bcf(y,x+1,z+1);
                if nb ==0
                   bctype = bct(y,x+1,z+1);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,10) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,10) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,11) = -k*G(j)/(4*dz*dx);
                nb = bcf(y,x-1,z+1);
                if nb ==0
                   bctype = bct(y,x-1,z+1);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1) ;
                      Mz(count,11) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,11) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,12) = -k*G(j)/(4*dz*dx);
                nb = bcf(y,x+1,z-1);
                if nb ==0
                   bctype = bct(y,x+1,z-1);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,12) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,12) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,13) = k*G(j)/(4*dz*dx);
                nb = bcf(y,x-1,z-1);
                if nb ==0
                   bctype = bct(y,x-1,z-1);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,13) = 0;
                   end
                elseif sum(nb == xtype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,13) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                
                Mz(count,14) = k*G(j)/(4*dz*dy);
                nb = bcf(y+1,x,z+1);
                if nb ==0
                   bctype = bct(y+1,x,z+1);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,14) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,14) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,15) = -k*G(j)/(4*dz*dy);
                nb = bcf(y-1,x,z+1);
                if nb ==0
                   bctype = bct(y-1,x,z+1);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,15) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,15) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,16) = -k*G(j)/(4*dz*dy);
                nb = bcf(y+1,x,z-1);
                if nb ==0
                   bctype = bct(y+1,x,z-1);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,16) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,16) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,17) = k*G(j)/(4*dz*dy);
                nb = bcf(y-1,x,z-1);
                if nb ==0
                   bctype = bct(y-1,x,z-1);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1) ;
                      Mz(count,17) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,17) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,18) = k*Gz(j)/(2*dy);
                nb = bcf(y+1,x,z);
                if nb ==0
                   bctype = bct(y+1,x,z);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,18) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,18) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Mz(count,19) = -k*Gz(j)/(2*dy);
                nb = bcf(y-1,x,z);
                if nb ==0
                   bctype = bct(y-1,x,z);
                   if bctype == 1 %dirchlett
                      Uz_bc(count,1) = Uz_bc(count,1);
                      Mz(count,19) = 0;
                   end
                elseif sum(nb == ytype) == 1
                    Uz_bc(count,1) = Uz_bc(count,1)+0;
                    Mz(count,19) = 0;
                else
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                end
                
                Fz_v(count,1) = Fz(j);
                F_list(count,1) = j;
                
            elseif f ~= 5 && f~= 0 && ft == 1
                count = count +1;
                M_loc(j) = count;    
                Ux_bc(count,1) = 0;
                Uy_bc(count,1) = 0;
                Uz_bc(count,1) = 0;
                
                % ---- Mx Build ------
                if sum(f == xtype) == 1
                    Mx(count,1) = 1;
                    Ux_bc(count,1) = 0;
                else
                    Mx(count,1) = -2*G(j)*((1/dx^2)+(1/dy^2)+(1/dz^2)+(k/dx^2));   
                    
                    Mx(count,2) = (Gx(j)/(2*dx))*(1+k)+(G(j)/dx^2)*(1+k);
                    nb = bcf(y,x+1,z);
                    if nb ==0
                        bctype  = bct(y,x+1,z);
                        if bctype == 1
                            Ux_bc(count,1) = Ux_bc(count,1);
                            Mx(count,2) = 0;
                        end
                    elseif sum(nb == xtype) == 1
                        Ux_bc(count,1) = Ux_bc(count,1)+0;
                        Mx(count,2) = 0;
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end                

                    Mx(count,3) = -(Gx(j)/(2*dx))*(1+k)+(G(j)/dx^2)*(1+k);
                    nb = bcf(y,x-1,z);
                    if nb ==0
                        bctype  = bct(y,x-1,z);
                        if bctype == 1
                            Ux_bc(count,1) = Ux_bc(count,1);
                            Mx(count,3) = 0;
                        end
                    elseif sum(nb == xtype) == 1
                        Ux_bc(count,1) = Ux_bc(count,1)+0;
                        Mx(count,3) = 0;
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end                
              
                    if sum(f == ytype) == 1
                        try
                            nb = bcf(y+1,x,z);
                            if nb ~= 0
                                Mx(count,4) =  (Gy(j)/(2*dy)) + G(j)/dy^2;
                                Mx(count,8) = k*Gx(j)/(2*dy);
                            end
                        catch
                            Mx(count,4) = 0;
                            Mx(count,8) = 0;
                        end
                        try
                            nb = bcf(y-1,x,z);
                            if nb ~= 0
                                Mx(count,5) = (-Gy(j)/(2*dy)) + G(j)/dy^2;
                                Mx(count,9) = -k*Gx(j)/(2*dy);
                            end
                        catch
                            Mx(count,5) = 0;
                            Mx(count,9) = 0;
                        end            

                    else
                        Mx(count,4) = (Gy(j)/(2*dy)) + G(j)/dy^2;
                        nb = bcf(y+1,x,z);
                        if nb ==0
                            bctype  = bct(y+1,x,z);
                            if bctype == 1
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,4) = 0;
                            end
                        elseif sum(nb == xtype) == 1
                            Ux_bc(count,1) = Ux_bc(count,1)+0;
                            Mx(count,4) = 0;
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end   
                        
                        Mx(count,5) = (-Gy(j)/(2*dy)) + G(j)/dy^2;
                        nb = bcf(y-1,x,z);
                        if nb ==0
                           bctype = bct(y-1,x,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,5) = 0;
                           end
                        elseif sum(nb == xtype) == 1
                            Ux_bc(count,1) = Ux_bc(count,1)+0;
                            Mx(count,5) = 0;
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end      

                        Mx(count,8) = k*Gx(j)/(2*dy);
                        nb = bcf(y+1,x,z);
                        if nb ==0
                           bctype = bct(y+1,x,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,8) = 0;
                           end
                        elseif sum(nb == ytype) == 1
                            Ux_bc(count,1) = Ux_bc(count,1)+0;
                            Mx(count,8) = 0;
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                

                        Mx(count,9) = -k*Gx(j)/(2*dy);
                        nb = bcf(y-1,x,z);
                        if nb ==0
                           bctype = bct(y-1,x,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,9) = 0;
                           end
                        elseif sum(nb == ytype) == 1
                            Ux_bc(count,1) = Ux_bc(count,1)+0;
                            Mx(count,9) = 0;
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                        

                    end

                    if sum(f==ztype) == 1
                        try
                            nb = bcf(y,x,z+1);
                            if nb ~= 0
                                Mx(count,6) = (Gz(j)/(2*dz)) + G(j)/dz^2;
                                Mx(count,14) = Gx(j)*k/(2*dz);
                            end
                        catch
                            Mx(count,[6 14]) = 0;
                        end
                        try
                            nb = bcf(y,x,z-1);
                            if nb ~= 0
                                Mx(count,7) = (-Gz(j)/(2*dz)) + G(j)/dz^2;
                                Mx(count,15) = -Gx(j)*k/(2*dz);
                            end
                        catch
                            Mx(count,[7 15]) = 0;
                        end     

                    else
                        Mx(count,6) = (Gz(j)/(2*dz)) + G(j)/dz^2;
                        nb = bcf(y,x,z+1);                  
                        if nb ==0
                           bctype = bct(y,x,z+1);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,6) = 0;
                           end
                        elseif sum(nb == xtype) == 1
                            Ux_bc(count,1) = Ux_bc(count,1)+0;
                            Mx(count,6) = 0;
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end   

                        Mx(count,7) = (-Gz(j)/(2*dz)) + G(j)/dz^2;
                        nb = bcf(y,x,z-1);
                        if nb ==0
                           bctype = bct(y,x,z-1);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,7) = 0;
                           end
                        elseif sum(nb == xtype) == 1
                            Ux_bc(count,1) = Ux_bc(count,1)+0;
                            Mx(count,7) = 0;
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end  

                        Mx(count,14) = Gx(j)*k/(2*dz);
                        nb = bcf(y,x,z+1);
                        if nb ==0
                           bctype = bct(y,x,z+1);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,14) = 0;
                           end
                        elseif sum(nb == ztype) == 1
                            Ux_bc(count,1) = Ux_bc(count,1)+0;
                            Mx(count,14) = 0;
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                

                        Mx(count,15) = -Gx(j)*k/(2*dz);
                        nb = bcf(y,x,z-1);
                        if nb ==0
                           bctype = bct(y,x,z-1);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,15) = 0;
                           end
                        elseif sum(nb == ztype) == 1
                            Ux_bc(count,1) = Ux_bc(count,1)+0;
                            Mx(count,15) = 0;
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end      
                    end

                    try
                        nb = bcf(y+1,x+1,z);
                        if nb ~= 0
                            Mx(count,10) = k*G(j)/(4*dx*dy);
                        end
                    catch
                        Mx(count,10) = 0;
                    end                             
       
                    try
                        nb = bcf(y-1,x+1,z);
                        if nb ~= 0
                            Mx(count,11) = -k*G(j)/(4*dx*dy);
                        end
                    catch
                        Mx(count,11) = 0;
                    end                             
                      
                    try
                        nb = bcf(y+1,x-1,z);
                        if nb ~= 0
                            Mx(count,12) = -k*G(j)/(4*dx*dy);
                        end
                    catch
                        Mx(count,12) = 0;
                    end                             
       
                    try
                        nb = bcf(y-1,x-1,z);
                        if nb ~= 0
                            Mx(count,13) = k*G(j)/(4*dx*dy);
                        end
                    catch
                        Mx(count,13) = 0;
                    end            
                  
  
                    try
                        nb = bcf(y,x+1,z+1);
                        if nb ~= 0
                            Mx(count,16) = k*G(j)/(4*dx*dz);
                        end
                    catch
                        Mx(count,16) = 0;
                    end      
           
                    try
                        nb = bcf(y,x+1,z-1);
                        if nb ~= 0
                            Mx(count,17) = -k*G(j)/(4*dx*dz);
                        end
                    catch
                        Mx(count,17) = 0;
                    end        
                
                    try
                        nb = bcf(y,x-1,z+1);
                        if nb ~= 0
                            Mx(count,18) = -k*G(j)/(4*dx*dz);
                        end
                    catch
                        Mx(count,18) = 0;
                    end      
           
                    try
                        nb = bcf(y,x-1,z-1);
                        if nb ~= 0
                            Mx(count,19) = k*G(j)/(4*dx*dz);
                        end
                    catch
                        Mx(count,19) = 0;
                    end                      
                end
                
                Fx_v(count,1) = Fx(j);
  
                % ---- My Build ------
                if sum(f == ytype) == 1                  
                    My(count,1) = 1;
                    Uy_bc(count,1) = Uy_bc(count,1) + 0;
                else                            
                    My(count,1) = -2*G(j)*((1/dx^2)+(1/dy^2)+(1/dz^2)+(k/dy^2));
                    My(count,2) = (Gy(j)/(2*dy))*(1+k)+(G(j)/dy^2)*(1+k);
                    nb = bcf(y+1,x,z);
                    if nb ==0
                       bctype = bct(y+1,x,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1) ;
                          My(count,2) = 0;
                       end
                    elseif sum(nb == ytype) == 1
                        Uy_bc(count,1) = Uy_bc(count,1)+0;
                        My(count,2) = 0;
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end
                
                    My(count,3) = -(Gy(j)/(2*dy))*(1+k)+(G(j)/dy^2)*(1+k);
                    nb = bcf(y-1,x,z);
                    if nb ==0
                       bctype = bct(y-1,x,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1) ;
                          My(count,3) = 0;
                       end
                    elseif sum(nb == ytype) == 1
                        Uy_bc(count,1) = Uy_bc(count,1)+0;
                        My(count,3) = 0;
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end 
                    
                    if sum(f == xtype) == 1
                        try
                            nb = bcf(y,x+1,z);
                            if nb ~= 0
                                My(count,4) =  (Gx(j)/(2*dx))+G(j)/dx^2;
                                My(count,8) =   k*Gy(j)/(2*dx);
                            end
                        catch
                            My(count,4) = 0;
                            My(count,8) = 0;
                        end
                        try
                            nb = bcf(y,x-1,z);
                            if nb ~= 0
                                My(count,5) = -(Gx(j)/(2*dx))+G(j)/dx^2;
                                My(count,9) =  -k*Gy(j)/(2*dx);
                            end
                        catch
                            My(count,5) = 0;
                            My(count,9) = 0;
                        end       
                        
                        
                    else
                        My(count,4) = (Gx(j)/(2*dx))+G(j)/dx^2;
                        nb = bcf(y,x+1,z);
                        if nb ==0
                           bctype = bct(y,x+1,z);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,4) = 0;
                           end
                        elseif sum(nb == ytype) == 1
                            Uy_bc(count,1) = Uy_bc(count,1)+0;
                            My(count,4) = 0;
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end                     

                        My(count,5) = -(Gx(j)/(2*dx))+G(j)/dx^2;
                        nb = bcf(y,x-1,z);
                        if nb ==0
                           bctype = bct(y,x-1,z);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,5) = 0;
                           end
                        elseif sum(nb == ytype) == 1
                            Uy_bc(count,1) = Uy_bc(count,1)+0;
                            My(count,5) = 0;
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end     
                        
                        My(count,8) = k*Gy(j)/(2*dx);
                        nb = bcf(y,x+1,z);
                        if nb ==0
                           bctype = bct(y,x+1,z);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,8) = 0;
                           end
                        elseif sum(nb == xtype) == 1
                            Uy_bc(count,1) = Uy_bc(count,1)+0;
                            My(count,8) = 0;
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end                     

                        My(count,9) = -k*Gy(j)/(2*dx);
                        nb = bcf(y,x-1,z);
                        if nb ==0
                           bctype = bct(y,x-1,z);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,9) = 0;
                           end
                        elseif sum(nb == xtype) == 1
                            Uy_bc(count,1) = Uy_bc(count,1)+0;
                            My(count,9) = 0;
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end    
                        
                    end          
                    
                    if sum(f == ztype) == 1                        
                       try
                            nb = bcf(y,x,z+1);
                            if nb ~= 0
                                My(count,6) =  (Gz(j)/(2*dz))+G(j)/dz^2;
                                My(count,14) =   k*Gy(j)/(2*dz);
                            end
                        catch
                            My(count,6) = 0;
                            My(count,14) = 0;
                        end
                        try
                            nb = bcf(y,x,z-1);
                            if nb ~= 0
                                My(count,7) = -(Gz(j)/(2*dz))+G(j)/dz^2;
                                My(count,15) =  -k*Gy(j)/(2*dz);
                            end
                        catch
                            My(count,7) = 0;
                            My(count,15) = 0;
                        end                         
                        
                    else
                        My(count,6) = (Gz(j)/(2*dz))+G(j)/dz^2;
                        nb = bcf(y,x,z+1);
                        if nb ==0
                           bctype = bct(y,x,z+1);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1) ;
                              My(count,6) = 0;
                           end
                        elseif sum(nb == ytype) == 1
                            Uy_bc(count,1) = Uy_bc(count,1)+0;
                            My(count,6) = 0;
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end                     

                        My(count,7) = -(Gz(j)/(2*dz))+G(j)/dz^2;
                        nb = bcf(y,x,z-1);
                        if nb ==0
                           bctype = bct(y,x,z-1);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,7) = 0;
                           end
                        elseif sum(nb == ytype) == 1
                            Uy_bc(count,1) = Uy_bc(count,1)+0;
                            My(count,7) = 0;
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end   
                        
                        My(count,14) = k*Gy(j)/(2*dz);
                        nb = bcf(y,x,z+1);
                        if nb ==0
                           bctype = bct(y,x,z+1);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,14) = 0;
                           end
                        elseif sum(nb == ztype) == 1
                            Uy_bc(count,1) = Uy_bc(count,1)+0;
                            My(count,14) = 0;
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end                     

                        My(count,15) = -k*Gy(j)/(2*dz);
                        nb = bcf(y,x,z-1);
                        if nb ==0
                           bctype = bct(y,x,z-1);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,15) = 0;
                           end
                        elseif sum(nb == ztype) == 1
                            Uy_bc(count,1) = Uy_bc(count,1)+0;
                            My(count,15) = 0;
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end   
                    end                    
                    
                    try
                        nb = bcf(y+1,x+1,z);
                        if nb ~= 0
                            My(count,10) = k*G(j)/(4*dx*dy);
                        end
                    catch
                        My(count,10) = 0;
                    end   

                    try
                        nb = bcf(y+1,x-1,z);
                        if nb ~= 0
                            My(count,11) = -k*G(j)/(4*dx*dy);
                        end
                    catch
                        My(count,11) = 0;
                    end                   

                    try
                        nb = bcf(y-1,x+1,z);
                        if nb ~= 0
                            My(count,12) = -k*G(j)/(4*dx*dy);
                        end
                    catch
                        My(count,12) = 0;
                    end   

                    try
                        nb = bcf(y-1,x-1,z);
                        if nb ~= 0
                            My(count,13) = k*G(j)/(4*dx*dy);
                        end
                    catch
                        My(count,13) = 0;
                    end                  

                
                    try
                        nb = bcf(y+1,x,z+1);
                        if nb ~= 0
                            My(count,16) = k*G(j)/(4*dz*dy);
                        end
                    catch
                        My(count,16) = 0;
                    end   
                
                    try
                        nb = bcf(y+1,x,z-1);
                        if nb ~= 0
                            My(count,17) = -k*G(j)/(4*dz*dy);
                        end
                    catch
                        My(count,17) = 0;
                    end                   

                    try
                        nb = bcf(y-1,x,z+1);
                        if nb ~= 0
                            My(count,18) = -k*G(j)/(4*dz*dy);
                        end
                    catch
                        My(count,18) = 0;
                    end   
                
                    try
                        nb = bcf(y-1,x,z-1);
                        if nb ~= 0
                            My(count,19) = k*G(j)/(4*dz*dy);
                        end
                    catch
                        My(count,19) = 0;
                    end     
                    
                end                              
                Fy_v(count,1) = Fy(j);
                
  
                % ---- Mz Build ------
                if sum(f == ztype) == 1
                    Mz(count,1) = 1;
                    Uz_bc(count,1) = Uz_bc(count,1) + 0;
                else             
                    Mz(count,1) = -2*G(j)*((1/dx^2)+(1/dy^2)+(1/dz^2)+(k/dz^2));
                    
                    Mz(count,6) = (Gz(j)/(2*dz))*(1+k)+(G(j)/dz^2)*(1+k);
                    nb = bcf(y,x,z+1);
                    if nb ==0
                       bctype = bct(y,x,z+1);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,6) = 0;
                       end
                    elseif sum(nb == ztype) == 1
                        Uz_bc(count,1) = Uz_bc(count,1)+0;
                        Mz(count,6) = 0;
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end   

                    Mz(count,7) = -(Gz(j)/(2*dz))*(1+k)+(G(j)/dz^2)*(1+k);
                    nb = bcf(y,x,z-1);
                    if nb ==0
                       bctype = bct(y,x,z-1);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,7) = 0;
                       end
                    elseif sum(nb == ztype) == 1
                        Uz_bc(count,1) = Uz_bc(count,1)+0;
                        Mz(count,7) = 0;
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end   
                    
                    if sum(f == xtype) == 1
                        try
                            nb = bcf(y,x+1,z);
                            if nb ~= 0
                                Mz(count,2) =  (G(j)/(dx^2))+Gx(j)/(2*dx);
                                Mz(count,8) = k*Gz(j)/(2*dx);  
                            end
                        catch
                            Mz(count,2) = 0;
                            Mz(count,8) = 0;
                        end
                        
                        try
                            nb = bcf(y,x-1,z);
                            if nb ~= 0
                                Mz(count,3) =  (G(j)/(dx^2))-Gx(j)/(2*dx);
                                Mz(count,9) = -k*Gz(j)/(2*dx);  
                            end
                        catch
                            Mz(count,3) = 0;
                            Mz(count,9) = 0;
                        end                        
                        
                        
                    else
                        Mz(count,2) = (G(j)/(dx^2))+Gx(j)/(2*dx);
                        nb = bcf(y,x+1,z);
                        if nb ==0
                           bctype = bct(y,x+1,z);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,2) = 0;
                           end
                        elseif sum(nb == ztype) == 1
                            Uz_bc(count,1) = Uz_bc(count,1)+0;
                            Mz(count,2) = 0;
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end

                        Mz(count,3) = (G(j)/(dx^2))-Gx(j)/(2*dx);
                        nb = bcf(y,x-1,z);
                        if nb ==0
                           bctype = bct(y,x-1,z);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,3) = 0;
                           end
                        elseif sum(nb == ztype) == 1
                            Uz_bc(count,1) = Uz_bc(count,1)+0;
                            Mz(count,3) = 0;
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end  
                        
                        Mz(count,8) = k*Gz(j)/(2*dx);
                        nb = bcf(y,x+1,z);
                        if nb ==0
                           bctype = bct(y,x+1,z);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,8) = 0;
                           end
                        elseif sum(nb == xtype) == 1
                            Uz_bc(count,1) = Uz_bc(count,1)+0;
                            Mz(count,8) = 0;
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end

                        Mz(count,9) = -k*Gz(j)/(2*dx);
                        nb = bcf(y,x-1,z);
                        if nb ==0
                           bctype = bct(y,x-1,z);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,9) = 0;
                           end
                        elseif sum(nb == xtype) == 1
                            Uz_bc(count,1) = Uz_bc(count,1)+0;
                            Mz(count,9) = 0;
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end
                        
                    end

                    if sum(f == ytype) == 1
                        try
                            nb = bcf(y+1,x,z);
                            if nb ~= 0
                                Mz(count,4) =  (G(j)/(dy^2))+Gy(j)/(2*dy);
                                Mz(count,18) = k*Gz(j)/(2*dy);
                            end
                        catch
                            Mz(count,4) = 0;
                            Mz(count,18) = 0;
                        end                       
                        
                        try
                            nb = bcf(y-1,x,z);
                            if nb ~= 0
                                Mz(count,5) =  (G(j)/(dy^2))-Gy(j)/(2*dy);
                                Mz(count,19) = -k*Gz(j)/(2*dy);
                            end
                        catch
                            Mz(count,5) = 0;
                            Mz(count,19) = 0;
                        end                             
                    else
                        Mz(count,4) = (G(j)/(dy^2))+Gy(j)/(2*dy);
                        nb = bcf(y+1,x,z);
                        if nb ==0
                           bctype = bct(y+1,x,z);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,4) = 0;
                           end
                        elseif sum(nb == ztype) == 1
                            Uz_bc(count,1) = Uz_bc(count,1)+0;
                            Mz(count,4) = 0;
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end                

                        Mz(count,5) = (G(j)/(dy^2))-Gy(j)/(2*dy);
                        nb = bcf(y-1,x,z);
                        if nb ==0
                           bctype = bct(y-1,x,z);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,5) = 0;
                           end
                        elseif sum(nb == ztype) == 1
                            Uz_bc(count,1) = Uz_bc(count,1)+0;
                            Mz(count,5) = 0;
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end                           

                        Mz(count,18) = k*Gz(j)/(2*dy);
                        nb = bcf(y+1,x,z);
                        if nb ==0
                           bctype = bct(y+1,x,z);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,18) = 0;
                           end
                        elseif sum(nb == ytype) == 1
                            Uz_bc(count,1) = Uz_bc(count,1)+0;
                            Mz(count,18) = 0;
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end

                        Mz(count,19) = -k*Gz(j)/(2*dy);
                        nb = bcf(y-1,x,z);
                        if nb ==0
                           bctype = bct(y-1,x,z);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,19) = 0;
                           end
                        elseif sum(nb == ytype) == 1
                            Uz_bc(count,1) = Uz_bc(count,1)+0;
                            Mz(count,19) = 0;
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end
                        
                    end

                    try
                        nb = bcf(y,x+1,z+1);
                        if nb ~= 0
                            Mz(count,10) = k*G(j)/(4*dx*dz);
                        end
                    catch
                        Mz(count,10) = 0;
                    end  
  
                    try
                        nb = bcf(y,x-1,z+1);
                        if nb ~= 0
                            Mz(count,11) = -k*G(j)/(4*dx*dz);
                        end
                    catch
                        Mz(count,11) = 0;
                    end                  

                    try
                        nb = bcf(y,x+1,z-1);
                        if nb ~= 0
                            Mz(count,12) = -k*G(j)/(4*dx*dz);
                        end
                    catch
                        Mz(count,12) = 0;
                    end  
  
                    try
                        nb = bcf(y,x-1,z-1);
                        if nb ~= 0
                            Mz(count,13) = k*G(j)/(4*dx*dz);
                        end
                    catch
                        Mz(count,13) = 0;
                    end      

                    try
                        nb = bcf(y+1,x,z+1);
                        if nb ~= 0
                            Mz(count,14) = k*G(j)/(4*dy*dz);
                        end
                    catch
                        Mz(count,14) = 0;
                    end  

                    try
                        nb = bcf(y-1,x,z+1);
                        if nb ~= 0
                            Mz(count,15) = -k*G(j)/(4*dy*dz);
                        end
                    catch
                        Mz(count,15) = 0;
                    end  
                
                    try
                        nb = bcf(y+1,x,z-1);
                        if nb ~= 0
                            Mz(count,16) = -k*G(j)/(4*dy*dz);
                        end
                    catch
                        Mz(count,16) = 0;
                    end  

                    try
                        nb = bcf(y-1,x,z-1);
                        if nb ~= 0
                            Mz(count,17) = k*G(j)/(4*dy*dz);
                        end
                    catch
                        Mz(count,17) = 0;
                    end  
                             

                end
                Fz_v(count,1) = Fz(j);
                F_list(count,1) = j;
                                
            elseif ft == 2  % No stress boundary conditions
                count = count +1;
                M_loc(j) = count;
                
                % -----  Mx Build
                Mx(count,1) = -2*G(j)*((1/dx^2)+(1/dy^2)+(1/dz^2)+(k/dx^2));
                
                if sum(f == fym) == 1   % Y - 1 missing
                    Mx(count,4) = 2*G(j)/dy^2;
                    nb = bcf(y+1,x,z);
                    if nb ~= 5
                       bctype = bct(y+1,x,z);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,4) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end 
                    
                    Mx(count,4) = 0; Mx(count,8) = 0; Mx(count,9) = 0;

                elseif sum(f==fyp) == 1 % Y + 1 missing
                    Mx(count,5) = 2*G(j)/dy^2;
                    nb = bcf(y-1,x,z);
                    if nb ~= 5
                       bctype = bct(y-1,x,z);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,5) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end  
                    
                    Mx(count,5) = 0; Mx(count,8) = 0; Mx(count,9) = 0;
                else
                    Mx(count,4) = (Gy(j)/(2*dy)) + G(j)/dy^2;
                    nb = bcf(y+1,x,z);
                    if nb ~= 5
                       bctype = bct(y+1,x,z);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,4) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end                

                    Mx(count,5) = (-Gy(j)/(2*dy)) + G(j)/dy^2;
                    nb = bcf(y-1,x,z);
                    if nb ~= 5
                       bctype = bct(y-1,x,z);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,5) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end  
                    
                    Mx(count,8) = k*Gx(j)/(2*dy);
                    nb = bcf(y+1,x,z);
                    if nb ~= 5
                       bctype = bct(y+1,x,z);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,8) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end                

                    Mx(count,9) = -k*Gx(j)/(2*dy);
                    nb = bcf(y-1,x,z);
                    if nb ~= 5
                       bctype = bct(y-1,x,z);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,9) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end  
                end
                
                if sum(f == fzm) == 1   % z - 1 missing
                    Mx(count,6) = 2*G(j)/dz^2;
                    nb = bcf(y,x,z+1);
                    if nb ~= 5
                       bctype = bct(y,x,z+1);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,6) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end 
                    
                    Mx(count,7) = 0; Mx(count,14) = 0; Mx(count,15) = 0;
                    
                elseif sum(f==fzp) == 1 % z + 1 missing
                    Mx(count,7) = 2*G(j)/dz^2;
                    nb = bcf(y,x,z-1);
                    if nb ~= 5
                       bctype = bct(y,x,z-1);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,7) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end   
                    
                    Mx(count,6) = 0; Mx(count,14) = 0; Mx(count,15) = 0;
                else
                    Mx(count,6) = (Gz(j)/(2*dz)) + G(j)/dz^2;
                    nb = bcf(y,x,z+1);
                    if nb ~= 5
                       bctype = bct(y,x,z+1);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,6) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end          

                    Mx(count,7) = (-Gz(j)/(2*dz)) + G(j)/dz^2;
                    nb = bcf(y,x,z-1);
                    if nb ~= 5
                       bctype = bct(y,x,z-1);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,7) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end   
                    
                    Mx(count,14) = Gx(j)*k/(2*dz);
                    nb = bcf(y,x,z+1);
                    if nb ~= 5
                       bctype = bct(y,x,z+1);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,14) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end                

                    Mx(count,15) = -Gx(j)*k/(2*dz);
                    nb = bcf(y,x,z-1);
                    if nb ~= 5
                       bctype = bct(y,x,z-1);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,15) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end    
                  
                end
                
                if sum(f == fx) == 0 % no x boundary conditions
                   if sum(f == fy) == 0 % no y boundary conditions
                        Mx(count,10) = k*G(j)/(4*dx*dy);
                        nb = bcf(y+1,x+1,z);
                        if nb ~= 5
                           bctype = bct(y+1,x+1,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,10) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                

                        Mx(count,11) = -k*G(j)/(4*dx*dy);
                        nb = bcf(y-1,x+1,z);
                        if nb ~= 5
                           bctype = bct(y-1,x+1,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,11) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end       

                        Mx(count,12) = -k*G(j)/(4*dx*dy);
                        nb = bcf(y+1,x-1,z);
                        if nb ~= 5
                           bctype = bct(y+1,x-1,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,12) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                

                        Mx(count,13) = k*G(j)/(4*dx*dy);
                        nb = bcf(y-1,x-1,z);
                        if nb ~= 5
                           bctype = bct(y-1,x,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,13) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end      
                   end
                   
                   if sum(f == fz) == 0 % no z boundary conditions
                        Mx(count,16) = k*G(j)/(4*dx*dz);
                        nb = bcf(y,x+1,z+1);
                        if nb ~= 5
                           bctype = bct(y,x+1,z+1);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,16) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                

                        Mx(count,17) = -k*G(j)/(4*dx*dz);
                        nb = bcf(y,x+1,z-1);
                        if nb ~= 5
                           bctype = bct(y,x+1,z-1);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,17) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end       

                        Mx(count,18) = -k*G(j)/(4*dx*dz);
                        nb = bcf(y,x-1,z+1);
                        if nb ~= 5
                           bctype = bct(y,x-1,z+1);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,18) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                

                        Mx(count,19) = k*G(j)/(4*dx*dz);
                        nb = bcf(y,x-1,z-1);
                        if nb ~= 5
                           bctype = bct(y-1,x,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,19) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end
                   end
                end  
                 
                Fx_v(count,1) = Fx(j);

                if sum(f == fxm) == 1   % x - 1 missing
                    Mx(count,2) = 2*(G(j)/dx^2)*(1+k);
                    nb = bcf(y,x+1,z);
                    if nb ~= 5
                       bctype = bct(y,x+1,z);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,2) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end
                    
                    Mx(count,3) = 0;
                elseif sum(f==fxp) == 1 % x + 1 missing
                    Mx(count,3) = 2*(G(j)/dx^2)*(1+k);
                    nb = bcf(y,x-1,z);
                    if nb ~= 5
                       bctype = bct(y,x-1,z);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,3) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end
                    
                    Mx(count,2) = 0;
                else
                    Mx(count,2) = (Gx(j)/(2*dx))*(1+k)+(G(j)/dx^2)*(1+k);
                    nb = bcf(y,x+1,z);
                    if nb ~= 5
                       bctype = bct(y,x+1,z);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,2) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end

                    Mx(count,3) = -(Gx(j)/(2*dx))*(1+k)+(G(j)/dx^2)*(1+k);
                    nb = bcf(y,x-1,z);
                    if nb ~= 5
                       bctype = bct(y,x-1,z);
                       if bctype == 1 %dirchlett
                          Ux_bc(count,1) = Ux_bc(count,1);
                          Mx(count,3) = 0;
                       end
                    else
                        Ux_bc(count,1) = Ux_bc(count,1) + 0;
                    end
                
                  
                end
             
                if sum(f == fx) == 0 % no x boundary conditions
                   if sum(f == fy) == 0 % no y boundary conditions
                        Mx(count,10) = k*G(j)/(4*dx*dy);
                        nb = bcf(y+1,x+1,z);
                        if nb ~= 5
                           bctype = bct(y+1,x+1,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,10) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                

                        Mx(count,11) = -k*G(j)/(4*dx*dy);
                        nb = bcf(y-1,x+1,z);
                        if nb ~= 5
                           bctype = bct(y-1,x+1,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,11) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end       

                        Mx(count,12) = -k*G(j)/(4*dx*dy);
                        nb = bcf(y+1,x-1,z);
                        if nb ~= 5
                           bctype = bct(y+1,x-1,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,12) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                

                        Mx(count,13) = k*G(j)/(4*dx*dy);
                        nb = bcf(y-1,x-1,z);
                        if nb ~= 5
                           bctype = bct(y-1,x,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,13) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end      
                   end
                   
                   if sum(f == fz) == 0 % no z boundary conditions
                        Mx(count,16) = k*G(j)/(4*dx*dz);
                        nb = bcf(y,x+1,z+1);
                        if nb ~= 5
                           bctype = bct(y,x+1,z+1);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,16) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                

                        Mx(count,17) = -k*G(j)/(4*dx*dz);
                        nb = bcf(y,x+1,z-1);
                        if nb ~= 5
                           bctype = bct(y,x+1,z-1);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,17) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end       

                        Mx(count,18) = -k*G(j)/(4*dx*dz);
                        nb = bcf(y,x-1,z+1);
                        if nb ~= 5
                           bctype = bct(y,x-1,z+1);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,18) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end                

                        Mx(count,19) = k*G(j)/(4*dx*dz);
                        nb = bcf(y,x-1,z-1);
                        if nb ~= 5
                           bctype = bct(y-1,x,z);
                           if bctype == 1 %dirchlett
                              Ux_bc(count,1) = Ux_bc(count,1);
                              Mx(count,19) = 0;
                           end
                        else
                            Ux_bc(count,1) = Ux_bc(count,1) + 0;
                        end
                   end
                end  
                 
                Fx_v(count,1) = Fx(j);

                % -----  My Build                
                My(count,1) = -2*G(j)*((1/dx^2)+(1/dy^2)+(1/dz^2)+(k/dy^2));
                
                if sum(f == fym) == 1   % Y - 1 missing
                    My(count,2) = (Gy(j)/(2*dy))*(1+k)+(G(j)/dy^2)*(1+k);
                    nb = bcf(y+1,x,z);
                    if nb ~= 5
                       bctype = bct(y+1,x,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,2) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end
                    
                    My(count,3) = 0;

                elseif sum(f==fyp) == 1 % Y + 1 missing
                    My(count,3) = -(Gy(j)/(2*dy))*(1+k)+(G(j)/dy^2)*(1+k);
                    nb = bcf(y-1,x,z);
                    if nb ~= 5
                       bctype = bct(y-1,x,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,3) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end  
                    
                    My(count,2) = 0;
                else
                    My(count,2) = (Gy(j)/(2*dy))*(1+k)+(G(j)/dy^2)*(1+k);
                    nb = bcf(y+1,x,z);
                    if nb ~= 5
                       bctype = bct(y+1,x,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,2) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end

                    My(count,3) = -(Gy(j)/(2*dy))*(1+k)+(G(j)/dy^2)*(1+k);
                    nb = bcf(y-1,x,z);
                    if nb ~= 5
                       bctype = bct(y-1,x,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,3) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end                                           

                end
                
                if sum(f == fzm) == 1   % z - 1 missing
                    My(count,6) = 2*G(j)/dz^2;
                    nb = bcf(y,x,z+1);
                    if nb ~= 5
                       bctype = bct(y,x,z+1);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,6) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end    
                    
                    My(count,7) = 0; My(count,14) = 0; My(count,15) = 0;
                    
                elseif sum(f==fzp) == 1 % z + 1 missing
                    My(count,7) = 2*G(j)/dz^2;
                    nb = bcf(y,x,z-1);
                    if nb ~= 5
                       bctype = bct(y,x,z-1);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,7) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end   
                    
                    My(count,6) = 0; My(count,14) = 0; My(count,15) = 0;
                else
                    My(count,6) = (Gz(j)/(2*dz))+G(j)/dz^2;
                    nb = bcf(y,x,z+1);
                    if nb ~= 5
                       bctype = bct(y,x,z+1);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,6) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end                     

                    My(count,7) = -(Gz(j)/(2*dz))+G(j)/dz^2;
                    nb = bcf(y,x,z-1);
                    if nb ~= 5
                       bctype = bct(y,x,z-1);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,7) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end   
                  
                    My(count,14) = k*Gy(j)/(2*dz);
                    nb = bcf(y,x,z+1);
                    if nb ~= 5
                       bctype = bct(y,x,z+1);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,14) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end                     

                    My(count,15) = -k*Gy(j)/(2*dz);
                    nb = bcf(y,x,z-1);
                    if nb ~= 5
                       bctype = bct(y,x,z-1);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,15) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end   
                end
                
                if sum(f == fxm) == 1   % x - 1 missing
                    My(count,4) = 2*G(j)/dx^2;
                    nb = bcf(y,x+1,z);
                    if nb ~= 5
                       bctype = bct(y,x+1,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,4) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end   
                    
                    My(count,5) = 0; My(count,8) = 0; My(count,9) = 0;

                elseif sum(f==fxp) == 1 % x + 1 missing
                    My(count,5) = 2*G(j)/dx^2;
                    nb = bcf(y,x-1,z);
                    if nb ~= 5
                       bctype = bct(y,x-1,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,5) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end    
                    
                    My(count,4) = 0; My(count,8) = 0; My(count,9) = 0;
                else
                    My(count,4) = (Gx(j)/(2*dx))+G(j)/dx^2;
                    nb = bcf(y,x+1,z);
                    if nb ~= 5
                       bctype = bct(y,x+1,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,4) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end                     

                    My(count,5) = -(Gx(j)/(2*dx))+G(j)/dx^2;
                    nb = bcf(y,x-1,z);
                    if nb ~= 5
                       bctype = bct(y,x-1,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,5) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end    

                    My(count,8) = k*Gy(j)/(2*dx);
                    nb = bcf(y,x+1,z);
                    if nb ~= 5
                       bctype = bct(y,x+1,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,8) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end                     

                    My(count,9) = -k*Gy(j)/(2*dx);
                    nb = bcf(y,x-1,z);
                    if nb ~= 5
                       bctype = bct(y,x-1,z);
                       if bctype == 1 %dirchlett
                          Uy_bc(count,1) = Uy_bc(count,1);
                          My(count,9) = 0;
                       end
                    else
                        Uy_bc(count,1) = Uy_bc(count,1) + 0;
                    end   

                  
                end
             
                if sum(f == fy) == 0 % no y boundary conditions
                   if sum(f == fx) == 0 % no x boundary conditions
                        My(count,10) = k*G(j)/(4*dx*dy);
                        nb = bcf(y+1,x+1,z);
                        if nb ~= 5
                           bctype = bct(y+1,x+1,z);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,10) = 0;
                           end
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end                     

                        My(count,11) = -k*G(j)/(4*dx*dy);
                        nb = bcf(y+1,x-1,z);
                        if nb ~= 5
                           bctype = bct(y+1,x-1,z);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,11) = 0;
                           end
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end  

                        My(count,12) = -k*G(j)/(4*dx*dy);
                        nb = bcf(y-1,x+1,z);
                        if nb ~= 5
                           bctype = bct(y-1,x+1,z);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,12) = 0;
                           end
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end                     

                        My(count,13) = k*G(j)/(4*dx*dy);
                        nb = bcf(y-1,x-1,z);
                        if nb ~= 5
                           bctype = bct(y-1,x-1,z);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,13) = 0;
                           end
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end  
                   end
                   
                   if sum(f == fz) == 0 % no z boundary conditions                
                        My(count,16) = k*G(j)/(4*dz*dy);
                        nb = bcf(y+1,x,z+1);
                        if nb ~= 5
                           bctype = bct(y+1,x,z+1);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,16) = 0;
                           end
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end                     

                        My(count,17) = -k*G(j)/(4*dz*dy);
                        nb = bcf(y+1,x,z-1);
                        if nb ~= 5
                           bctype = bct(y+1,x,z-1);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,17) = 0;
                           end
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end   

                        My(count,18) = -k*G(j)/(4*dz*dy);
                        nb = bcf(y-1,x,z+1);
                        if nb ~= 5
                           bctype = bct(y-1,x,z+1);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,18) = 0;
                           end
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end                     

                        My(count,19) = k*G(j)/(4*dz*dy);
                        nb = bcf(y-1,x,z-1);
                        if nb ~= 5
                           bctype = bct(y-1,x,z-1);
                           if bctype == 1 %dirchlett
                              Uy_bc(count,1) = Uy_bc(count,1);
                              My(count,19) = 0;
                           end
                        else
                            Uy_bc(count,1) = Uy_bc(count,1) + 0;
                        end
                   end
                end  
                 
                Fy_v(count,1) = Fy(j);
                
                % -----  Mz Build     
                Mz(count,1) = -2*G(j)*((1/dx^2)+(1/dy^2)+(1/dz^2)+(k/dz^2));
                
                if sum(f == fym) == 1   % Y - 1 missing
                    Mz(count,4) = 2*(G(j)/(dy^2));
                    nb = bcf(y+1,x,z);
                    if nb ~= 5
                       bctype = bct(y+1,x,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,4) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end    
                    
                    Mz(count,5) = 0; Mz(count,18) = 0; Mz(count,19) = 0;

                elseif sum(f==fyp) == 1 % Y + 1 missing
                    Mz(count,5) = 2*(G(j)/(dy^2));
                    nb = bcf(y-1,x,z);
                    if nb ~= 5
                       bctype = bct(y-1,x,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,5) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end 
                    
                    Mz(count,4) = 0; Mz(count,18) = 0; Mz(count,19) = 0;

                else
                    Mz(count,4) = (G(j)/(dy^2))+Gy(j)/(2*dy);
                    nb = bcf(y+1,x,z);
                    if nb ~= 5
                       bctype = bct(y+1,x,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,4) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end                

                    Mz(count,5) = (G(j)/(dy^2))-Gy(j)/(2*dy);
                    nb = bcf(y-1,x,z);
                    if nb ~= 5
                       bctype = bct(y-1,x,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,5) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end    
                    
                    Mz(count,18) = k*Gz(j)/(2*dy);
                    nb = bcf(y+1,x,z);
                    if nb ~= 5
                       bctype = bct(y+1,x,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,18) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end

                    Mz(count,19) = -k*Gz(j)/(2*dy);
                    nb = bcf(y-1,x,z);
                    if nb ~= 5
                       bctype = bct(y-1,x,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,19) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end


                end
                
                if sum(f == fzm) == 1   % z - 1 missing
                    Mz(count,6) = 2*(G(j)/dz^2)*(1+k);
                    nb = bcf(y,x,z+1);
                    if nb ~= 5
                       bctype = bct(y,x,z+1);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,6) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end   
                    
                    Mz(count,7) = 0;
                    
                elseif sum(f==fzp) == 1 % z + 1 missing
                    Mz(count,7) = 2*(G(j)/dz^2)*(1+k);
                    nb = bcf(y,x,z-1);
                    if nb ~= 5
                       bctype = bct(y,x,z-1);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,7) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end  
                    
                    Mz(count,6) = 0;
                else
                    Mz(count,6) = (Gz(j)/(2*dz))*(1+k)+(G(j)/dz^2)*(1+k);
                    nb = bcf(y,x,z+1);
                    if nb ~= 5
                       bctype = bct(y,x,z+1);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,6) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end   

                    Mz(count,7) = -(Gz(j)/(2*dz))*(1+k)+(G(j)/dz^2)*(1+k);
                    nb = bcf(y,x,z-1);
                    if nb ~= 5
                       bctype = bct(y,x,z-1);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,7) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end   
                  

                end
                
                if sum(f == fxm) == 1   % x - 1 missing
                    Mz(count,2) = (G(j)/(dx^2))+Gx(j)/(2*dx);
                    nb = bcf(y,x+1,z);
                    if nb ~= 5
                       bctype = bct(y,x+1,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,2) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end
                    
                    Mz(count,3) = 0; Mz(count,8) = 0; Mz(count,9) = 0;

                elseif sum(f==fxp) == 1 % x + 1 missing
                    Mz(count,3) = (G(j)/(dx^2))-Gx(j)/(2*dx);
                    nb = bcf(y,x-1,z);
                    if nb ~= 5
                       bctype = bct(y,x-1,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,3) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end  
                    
                    Mz(count,2) = 0; Mz(count,8) = 0; Mz(count, 9) = 0;
                else
                    Mz(count,2) = (G(j)/(dx^2))+Gx(j)/(2*dx);
                    nb = bcf(y,x+1,z);
                    if nb ~= 5
                       bctype = bct(y,x+1,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,2) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end

                    Mz(count,3) = (G(j)/(dx^2))-Gx(j)/(2*dx);
                    nb = bcf(y,x-1,z);
                    if nb ~= 5
                       bctype = bct(y,x-1,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,3) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end   
                    
                    Mz(count,8) = k*Gz(j)/(2*dx);
                    nb = bcf(y,x+1,z);
                    if nb ~= 5
                       bctype = bct(y,x+1,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,8) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end

                    Mz(count,9) = -k*Gz(j)/(2*dx);
                    nb = bcf(y,x-1,z);
                    if nb ~= 5
                       bctype = bct(y,x-1,z);
                       if bctype == 1 %dirchlett
                          Uz_bc(count,1) = Uz_bc(count,1);
                          Mz(count,9) = 0;
                       end
                    else
                        Uz_bc(count,1) = Uz_bc(count,1) + 0;
                    end
                end
                
             
                if sum(f == fz) == 0 % no y boundary conditions
                   if sum(f == fx) == 0 % no x boundary conditions
                        Mz(count,10) = k*G(j)/(4*dz*dx);
                        nb = bcf(y,x+1,z+1);
                        if nb ~= 5
                           bctype = bct(y,x+1,z+1);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,10) = 0;
                           end
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end

                        Mz(count,11) = -k*G(j)/(4*dz*dx);
                        nb = bcf(y,x-1,z+1);
                        if nb ~= 5
                           bctype = bct(y,x-1,z+1);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,11) = 0;
                           end
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end

                        Mz(count,12) = -k*G(j)/(4*dz*dx);
                        nb = bcf(y,x+1,z-1);
                        if nb ~= 5
                           bctype = bct(y,x+1,z-1);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,12) = 0;
                           end
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end

                        Mz(count,13) = k*G(j)/(4*dz*dx);
                        nb = bcf(y,x-1,z-1);
                        if nb ~= 5
                           bctype = bct(y,x-1,z-1);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,13) = 0;
                           end
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end
                       
                   end
                   
                   if sum(f == fy) == 0 % no z boundary conditions                
                        Mz(count,14) = k*G(j)/(4*dz*dy);
                        nb = bcf(y+1,x,z+1);
                        if nb ~= 5
                           bctype = bct(y+1,x,z+1);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,14) = 0;
                           end
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end

                        Mz(count,15) = -k*G(j)/(4*dz*dy);
                        nb = bcf(y-1,x,z+1);
                        if nb ~= 5
                           bctype = bct(y-1,x,z+1);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,15) = 0;
                           end
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end

                        Mz(count,16) = -k*G(j)/(4*dz*dy);
                        nb = bcf(y+1,x,z-1);
                        if nb ~= 5
                           bctype = bct(y+1,x,z-1);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,16) = 0;
                           end
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end

                        Mz(count,17) = k*G(j)/(4*dz*dy);
                        nb = bcf(y-1,x,z-1);
                        if nb ~= 5
                           bctype = bct(y-1,x,z-1);
                           if bctype == 1 %dirchlett
                              Uz_bc(count,1) = Uz_bc(count,1);
                              Mz(count,17) = 0;
                           end
                        else
                            Uz_bc(count,1) = Uz_bc(count,1) + 0;
                        end  
                   end
                end  
                 
                Fz_v(count,1) = Fz(j);                              
                
            end
        end
    end
end




Fv = [Fx_v; Fy_v; Fz_v];

b = sum(Mx(:)~=0) + sum(My(:) ~= 0) + sum(Mz(:)~=0);
M_dx = zeros(sum(Mx(:)~=0),3);
M_dy = zeros(sum(My(:)~=0),3);
M_dz = zeros(sum(Mz(:)~=0),3);


M_loc_m = zeros(count,3);
j = 0; i = 0; kk = 0;
for z = 1:sz
    for y = 1:sy
        for x = 1:sx 

            if M_loc(y,x,z) ~= 0
               eq_line = M_loc(y,x,z);
               M_loc_m(eq_line,1) = y; M_loc_m(eq_line,2) = x; M_loc_m(eq_line,3) = z;
               
               % ---------   Mx   Build  ----------
               if Mx(eq_line,1) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc;
                   M_dx(j,3) = Mx(eq_line,1);                   
               end
               
               if Mx(eq_line,2) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x+1,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc;
                   M_dx(j,3) = Mx(eq_line,2);                   
               end
               
               if Mx(eq_line,3) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x-1,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc;
                   M_dx(j,3) = Mx(eq_line,3);                   
               end

               if Mx(eq_line,4) ~= 0
                   j = j + 1;
                   loc = M_loc(y+1,x,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc;
                   M_dx(j,3) = Mx(eq_line,4);                   
               end
               
               if Mx(eq_line,5) ~= 0
                   j = j + 1;
                   loc = M_loc(y-1,x,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc;
                   M_dx(j,3) = Mx(eq_line,5);                   
               end
               
               if Mx(eq_line,6) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x,z+1);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc;
                   M_dx(j,3) = Mx(eq_line,6);                   
               end
               
               if Mx(eq_line,7) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x,z-1);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc;
                   M_dx(j,3) = Mx(eq_line,7);                   
               end
       
               if Mx(eq_line,8) ~= 0
                   j = j + 1;
                   loc = M_loc(y+1,x,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + count;
                   M_dx(j,3) = Mx(eq_line,8);                   
               end
               
               if Mx(eq_line,9) ~= 0
                   j = j + 1;
                   loc = M_loc(y-1,x,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + count;
                   M_dx(j,3) = Mx(eq_line,9);                   
               end
                    
               if Mx(eq_line,10) ~= 0
                   j = j + 1;
                   loc = M_loc(y+1,x+1,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + count;
                   M_dx(j,3) = Mx(eq_line,10);                   
               end
               
               if Mx(eq_line,11) ~= 0
                   j = j + 1;
                   loc = M_loc(y-1,x+1,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + count;
                   M_dx(j,3) = Mx(eq_line,11);                   
               end   
               
               if Mx(eq_line,12) ~= 0
                   j = j + 1;
                   loc = M_loc(y+1,x-1,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + count;
                   M_dx(j,3) = Mx(eq_line,12);                   
               end
               
               if Mx(eq_line,13) ~= 0
                   j = j + 1;
                   loc = M_loc(y-1,x-1,z);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + count;
                   M_dx(j,3) = Mx(eq_line,13);                   
               end    

               if Mx(eq_line,14) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x,z+1);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + 2*count;
                   M_dx(j,3) = Mx(eq_line,14);                   
               end
               
               if Mx(eq_line,15) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x,z-1);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + 2*count;
                   M_dx(j,3) = Mx(eq_line,15);                   
               end                   
               
               if Mx(eq_line,16) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x+1,z+1);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + 2*count;
                   M_dx(j,3) = Mx(eq_line,16);                   
               end
               
               if Mx(eq_line,17) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x+1,z-1);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + 2*count;
                   M_dx(j,3) = Mx(eq_line,17);                   
               end
               
               if Mx(eq_line,18) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x-1,z+1);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + 2*count;
                   M_dx(j,3) = Mx(eq_line,18);                   
               end
               
               if Mx(eq_line,19) ~= 0
                   j = j + 1;
                   loc = M_loc(y,x-1,z-1);
                   M_dx(j,1) = eq_line;
                   M_dx(j,2) = loc + 2*count;
                   M_dx(j,3) = Mx(eq_line,19);                   
               end               
               
               % ---------   My   Build  ----------
               if My(eq_line,1) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y,x,z);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc + count;
                   M_dy(kk,3) = My(eq_line,1);                   
               end
               
               if My(eq_line,2) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y+1,x,z);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc + count;
                   M_dy(kk,3) = My(eq_line,2);                   
               end
               
               if My(eq_line,3) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y - 1,x,z);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc + count;
                   M_dy(kk,3) = My(eq_line,3);                   
               end

               if My(eq_line,4) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y,x+1,z);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc+count;
                   M_dy(kk,3) = My(eq_line,4);                   
               end
               
               if My(eq_line,5) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y,x-1,z);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc+ count;
                   M_dy(kk,3) = My(eq_line,5);                   
               end
               
               if My(eq_line,6) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y,x,z+1);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc+count;
                   M_dy(kk,3) = My(eq_line,6);                   
               end
               
               if My(eq_line,7) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y,x,z-1);
                   
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc+count;
                   M_dy(kk,3) = My(eq_line,7);                   
               end
       
               if My(eq_line,8) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y,x+1,z);
                                      if loc == 0
                       disp([y x z eq_line]);
                   end
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc;
                   M_dy(kk,3) = My(eq_line,8);                   
               end
               
               if My(eq_line,9) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y,x-1,z);
                                      if loc == 0
                       disp([y x z eq_line]);
                   end
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc;
                   M_dy(kk,3) = My(eq_line,9);                   
               end
                    
               if My(eq_line,10) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y+1,x+1,z);
                   if loc == 0
                       disp([y x z eq_line]);
                   end
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc;
                   M_dy(kk,3) = My(eq_line,10);                   
               end
               
               if My(eq_line,11) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y+1,x-1,z);
                                      if loc == 0
                       disp([y x z eq_line]);
                   end
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc;
                   M_dy(kk,3) = My(eq_line,11);                   
               end   
               
               if My(eq_line,12) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y-1,x+1,z);
                                      if loc == 0
                       disp([y x z eq_line]);
                   end
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc;
                   M_dy(kk,3) = My(eq_line,12);                   
               end
               
               if My(eq_line,13) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y-1,x-1,z);
                                      if loc == 0
                       disp([y x z eq_line]);
                   end
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc;
                   M_dy(kk,3) = My(eq_line,13);                   
               end    

               if My(eq_line,14) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y,x,z+1);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc + 2*count;
                   M_dy(kk,3) = My(eq_line,14);                   
               end
               
               if My(eq_line,15) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y,x,z-1);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc + 2*count;
                   M_dy(kk,3) = My(eq_line,15);                   
               end                   
               
               if My(eq_line,16) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y+1,x,z+1);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc + 2*count;
                   M_dy(kk,3) = My(eq_line,16);                   
               end
               
               if My(eq_line,17) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y+1,x,z-1);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc + 2*count;
                   M_dy(kk,3) = My(eq_line,17);                   
               end
               
               if My(eq_line,18) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y-1,x,z+1);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc + 2*count;
                   M_dy(kk,3) = My(eq_line,18);                   
               end
               
               if My(eq_line,19) ~= 0
                   kk= kk+ 1;
                   loc = M_loc(y-1,x,z-1);
                   M_dy(kk,1) = eq_line;
                   M_dy(kk,2) = loc + 2*count;
                   M_dy(kk,3) = My(eq_line,19);                   
               end               
               
               % ---------   Mz   Build  ----------
               if Mz(eq_line,1) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x,z);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc + 2*count;
                   M_dz(i,3) = Mz(eq_line,1);                   
               end
               
               if Mz(eq_line,2) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x+1,z);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc + 2*count;
                   M_dz(i,3) = Mz(eq_line,2);                   
               end
               
               if Mz(eq_line,3) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x-1,z);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc + 2*count;
                   M_dz(i,3) = Mz(eq_line,3);                   
               end

               if Mz(eq_line,4) ~= 0
                   i= i+ 1;
                   loc = M_loc(y+1,x,z);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc+2*count;
                   M_dz(i,3) = Mz(eq_line,4);                   
               end
               
               if Mz(eq_line,5) ~= 0
                   i= i+ 1;
                   loc = M_loc(y-1,x,z);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc+ 2*count;
                   M_dz(i,3) = Mz(eq_line,5);                   
               end
               
               if Mz(eq_line,6) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x,z+1);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc+2*count;
                   M_dz(i,3) = Mz(eq_line,6);                   
               end
               
               if Mz(eq_line,7) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x,z-1);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc+2*count;
                   M_dz(i,3) = Mz(eq_line,7);                   
               end
       
               if Mz(eq_line,8) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x+1,z);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc;
                   M_dz(i,3) = Mz(eq_line,8);                   
               end
               
               if Mz(eq_line,9) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x-1,z);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc;
                   M_dz(i,3) = Mz(eq_line,9);                   
               end
                    
               if Mz(eq_line,10) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x+1,z+1);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc;
                   M_dz(i,3) = Mz(eq_line,10);                   
               end
               
               if Mz(eq_line,11) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x-1,z+1);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc;
                   M_dz(i,3) = Mz(eq_line,11);                   
               end   
               
               if Mz(eq_line,12) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x+1,z-1);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc;
                   M_dz(i,3) = Mz(eq_line,12);                   
               end
               
               if Mz(eq_line,13) ~= 0
                   i= i+ 1;
                   loc = M_loc(y,x-1,z-1);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc;
                   M_dz(i,3) = Mz(eq_line,13);                   
               end    

               if Mz(eq_line,14) ~= 0
                   i= i+ 1;
                   loc = M_loc(y+1,x,z+1);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc + count;
                   M_dz(i,3) = Mz(eq_line,14);                   
               end
               
               if Mz(eq_line,15) ~= 0
                   i= i+ 1;
                   loc = M_loc(y-1,x,z+1);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc + count;
                   M_dz(i,3) = Mz(eq_line,15);                   
               end                   
               
               if Mz(eq_line,16) ~= 0
                   i= i+ 1;
                   loc = M_loc(y+1,x,z-1);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc + count;
                   M_dz(i,3) = Mz(eq_line,16);                   
               end
               
               if Mz(eq_line,17) ~= 0
                   i= i+ 1;
                   loc = M_loc(y-1,x,z-1);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc + count;
                   M_dz(i,3) = Mz(eq_line,17);                   
               end
               
               if Mz(eq_line,18) ~= 0
                   i= i+ 1;
                   loc = M_loc(y+1,x,z);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc + count;
                   M_dz(i,3) = Mz(eq_line,18);                   
               end
               
               if Mz(eq_line,19) ~= 0
                   i= i+ 1;
                   loc = M_loc(y-1,x,z);
                   M_dz(i,1) = eq_line;
                   M_dz(i,2) = loc + count;
                   M_dz(i,3) = Mz(eq_line,19);                   
               end               
            end
        end
    end 
end

M_dy(:,1) = M_dy(:,1)+count;
M_dz(:,1) = M_dz(:,1)+count*2;

totsnodes = sum(Mx(:)~=0)+sum(My(:)~=0)+sum(Mz(:)~=0);
M = zeros(totsnodes,3);
M(1:sum(Mx(:)~=0),:) = M_dx;
M(sum(Mx(:)~=0)+1:sum(Mx(:)~=0)+sum(My(:)~=0),:) = M_dy;
M(sum(Mx(:)~=0)+1+sum(My(:)~=0):sum(Mx(:)~=0)+sum(My(:)~=0)+sum(Mz(:)~=0),:) = M_dz;

Md = sparse(M(:,1),M(:,2),M(:,3),size(Fv,1),size(Fv,1));

end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end of file