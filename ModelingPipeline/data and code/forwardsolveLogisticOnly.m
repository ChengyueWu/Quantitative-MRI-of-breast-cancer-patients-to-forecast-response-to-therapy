%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Forward finite difference scheme
%Function receives variables pertaining to the diffusion and proliferation
%of tumor cells and returns the changes in the  output per the number of 
%time steps.

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [TCs] = forwardsolveLogisticOnly(TCs_in,dt_TC,k,carcap,bcf,sy_TC,sx_TC,sz_TC,Dtc,dx_TC,dy_TC,dz_TC,steps,mechstep)
% Inputs:
%      TCs_in   :     Tumor Cells In                  double  (sy,sx,sz)
%      dt_TC    :     time step                       double  (1x1)
%      k        :     growth rate                     double  (sy,sx,sz)
%      carcap   :     Carrying Capacity               double  (sy,sx,sz)
%      bcf      :     boundary conditions             double  (sy,sx,sz)
%      sx_TC    :     x matrix dimesion               double  (1x1)
%      sy_TC    :     y matrix dimension              double  (1x1)
%      sz_TC    :     z matrix dimension              double  (1x1)
%      Dtc      :     Tumor Cell Diffusion            double  (sy,sx,sz)
%      dx_TC    :     grid spacing in x               double  (1x1)
%      dy_TC    :     grid spacing in y               double  (1x1)
%      dz_TC    :     grid spacing in z               double  (1x1)
%      steps    :     step start                      double  (1x1)
%      mechsteps:     step end                        double  (1x1)
% Outputs:
%      TCs     :     Tumor Cells out                 double(sy,sx,sz)
% 
% The code for the boundaries is
% BCF = 5;   no boundaries
% BCF [1xx, 2xx, 3xx]:  a node is missing in the z-direction
% BCF [X1X X2X x3x]: a node is missing in the x direction
% BCF [xx1 xx2 xx3]: a node is missing in the y-direction.
% or
% 1 = missing a node above for z, ahead for x and under for y
% 2 = missing a node below for z, behind for x, and above for y
% 3 = missing in both directions
% else = any mixture of these three directions’ boundaries
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Holding matrices for the tumor cells and drug distribution
TCs = zeros(sy_TC,sx_TC,sz_TC);

for jj = steps:steps+mechstep
    %Calculate the logistic term and ratio of cells per carrying capacity
    growth_limit = (1-((TCs_in./carcap)));
    CC = TCs_in./carcap;
    %Nested loops to go through each voxel
    for z = 1:sz_TC
        for x = 1:sx_TC
            for y = 1:sy_TC
                f = bcf(y,x,z);
                if f ~= 0
                    %Calculate the right hand side of the
                    %reaction-diffusion equation
                    %Note that this code does not allow 2 negatives to
                    %results in positive growth
                    growth = k(y,x,z)*TCs_in(y,x,z)*growth_limit(y,x,z)*(k(y,x,z)>=0)...
                             + k(y,x,z)*TCs_in(y,x,z)*abs(growth_limit(y,x,z))*(k(y,x,z)<0);
                    if f == 5
                        %--- TC Deriv (1st deriv)
                        diff_1y = (1/(2*dy_TC))*(CC(y+1,x,z)-CC(y-1,x,z));
                        diff_1x = (1/(2*dx_TC))*(CC(y,x+1,z)-CC(y,x-1,z));
                        diff_1z = (1/(2*dz_TC))*(CC(y,x,z+1)-CC(y,x,z-1));
                        %--- Diffusion Constant Deriv (1st deriv)
                        diff_Dy = (1/(2*dy_TC))*(Dtc(y+1,x,z)-Dtc(y-1,x,z));
                        diff_Dx = (1/(2*dx_TC))*(Dtc(y,x+1,z)-Dtc(y,x-1,z));
                        diff_Dz = (1/(2*dz_TC))*(Dtc(y,x,z+1)-Dtc(y,x,z-1));
                        %--- TC Deriv (2nd deriv)
                        diff_y = (CC(y+1,x,z)-2*CC(y,x,z)+CC(y-1,x,z))/(dy_TC^2);
                        diff_x = (CC(y,x+1,z)-2*CC(y,x,z)+CC(y,x-1,z))/(dx_TC^2); 
                        diff_z = (CC(y,x,z+1)-2*CC(y,x,z)+CC(y,x,z-1))/(dz_TC^2); 

                    else
                        if f >= 300
                            f = f - 300;
                            diff_1z = 0;
                            diff_Dz = 0;
                            diff_z = 0;
                        elseif f >= 200
                            f = f - 200;
                            diff_z = (2*CC(y,x,z+1)-2*CC(y,x,z))/(dz_TC^2);
                            diff_1z = 0;
                            diff_Dz = 0;
                        elseif f >= 100
                            f = f - 100; 
                            diff_z = (2*CC(y,x,z-1)-2*CC(y,x,z))/(dz_TC^2);  
                            diff_1z = 0;
                            diff_Dz = 0;    
                        else
                            diff_z = (CC(y,x,z+1)-2*CC(y,x,z)+CC(y,x,z-1))/(dz_TC^2); 
                            diff_1z = (1/(2*dz_TC))*(CC(y,x,z+1)-CC(y,x,z-1));
                            diff_Dz = (1/(2*dz_TC))*(Dtc(y,x,z+1)-Dtc(y,x,z-1));
                        end

                        if f >= 30
                           f = f-30;
                            diff_1x = 0;
                            diff_Dx = 0;
                            diff_x = 0; 
                        elseif f >= 20
                            f = f-20;
                            diff_x = (2*CC(y,x+1,z)-2*CC(y,x,z))/(dy_TC^2); 
                            diff_Dx = 0;
                            diff_1x = 0;    
                        elseif f >= 10
                            f = f-10;
                            diff_x = (2*CC(y,x-1,z)-2*CC(y,x,z))/(dy_TC^2); 
                            diff_1x = 0;
                            diff_Dx = 0;
                        else
                            diff_x = (CC(y,x+1,z)-2*CC(y,x,z)+CC(y,x-1,z))/(dy_TC^2);
                            diff_1x = (1/(2*dx_TC))*(CC(y,x+1,z)-CC(y,x-1,z));
                            diff_Dx = (1/(2*dx_TC))*(Dtc(y,x+1,z)-Dtc(y,x-1,z));
                        end

                        if f == 3
                            diff_y = 0; 
                            diff_1y = 0;
                            diff_Dy = 0;
                        elseif f == 2
                            diff_y = (2*CC(y-1,x,z)-2*CC(y,x,z))/(dy_TC^2); 
                            diff_1y = 0;
                            diff_Dy = 0;
                        elseif f == 1 
                            diff_y = (2*CC(y+1,x,z)-2*CC(y,x,z))/(dy_TC^2);
                            diff_1y = 0;
                            diff_Dy = 0;
                        else
                            diff_y = (CC(y+1,x,z)-2*CC(y,x,z)+CC(y-1,x,z))/(dy_TC^2);
                            diff_Dy = (1/(2*dy_TC))*(Dtc(y+1,x,z)-Dtc(y-1,x,z));
                            diff_1y = (1/(2*dy_TC))*(CC(y+1,x,z)-CC(y-1,x,z));
                        end

                    end 
                    
                    %Calculate the next step
                    TCs(y,x,z) = TCs_in(y,x,z) + ...
                        dt_TC*(growth+Dtc(y,x,z)*(diff_y+diff_x+diff_z)+...
                        diff_Dx*diff_1x+diff_Dy*diff_1y+diff_Dz*diff_1z);      
                    if TCs(y,x,z) < 0 
                       TCs(y,x,z) = 0;
                    end
                else
                    TCs(y,x,z) = 0;
                end
            end %end of y dimesion loop
        end %end of x dimension loop
    end %end of z dimension loop
    
    %Update variable
    TCs_in = TCs;
    
end %end of time step loop

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end of file