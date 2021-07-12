%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Function to differentiating while also accounting for boundary conditions
%Using centered difference within the region and forward or backward
%difference for appropriate boundaries

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [Diffyx,Diffyy,Diffyz] = Diffy3D(dx,dy,dz,bcf,CC)
% Inputs:
%      dx_TC   :     grid spacing in x               double  (1x1)
%      dy_TC   :     grid spacing in y               double  (1x1)
%      dz_TC   :     grid spacing in z               double  (1x1)
%      bcf     :     boundary conditions             double  (sy,sx,sz)
%      CC      :     matrix to be differentiated     double  (sy,sx,sz)
% Outputs:
%      Diffyx  :     x derivative                    double  (sy,sx,sz)
%      Diffyy  :     y derivative                    double  (sy,sx,sz)
%      Diffyz  :     z derivative                    double  (sy,sx,sz)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[sy,sx,sz] = size(CC);
Diffyx = zeros(size(CC));
Diffyy = zeros(size(CC));
Diffyz = zeros(size(CC));

for z = 1:sz
    for x = 1:sx
        for y = 1:sy
            diff_1x = 0;
            diff_1y = 0;
            diff_1z = 0;
            f = bcf(y,x,z);
            if f ~= 0
                if f == 5
                    %TC Deriv (1st deriv)
                    diff_1y = (1/(2*dy))*(CC(y+1,x,z)-CC(y-1,x,z));
                    diff_1x = (1/(2*dx))*(CC(y,x+1,z)-CC(y,x-1,z));
                    diff_1z = (1/(2*dz))*(CC(y,x,z+1)-CC(y,x,z-1));
                else
                    if f >= 300
                        f = f - 300;
                        diff_1z = 0;
                    elseif f >= 200
                        f = f - 200;
                        diff_1z = 0;
                    elseif f >= 100
                        f = f - 100; 
                        diff_1z = 0;
                    else
                        diff_1z = (1/(2*dz))*(CC(y,x,z+1)-CC(y,x,z-1));
                    end

                    if f >= 30
                       f = f-30;
                        diff_1x = 0;
                    elseif f >= 20
                        f = f-20;
                        diff_1x = 0;    
                    elseif f >= 10
                        f = f-10;
                        diff_1x = 0;
                    else
                        diff_1x = (1/(2*dx))*(CC(y,x+1,z)-CC(y,x-1,z));
                    end

                    if f == 3
                        diff_1y = 0;
                    elseif f == 2
                        diff_1y = 0;
                    elseif f == 1 
                        diff_1y = 0;
                    else
                            diff_1y = (1/(2*dy))*(CC(y+1,x,z)-CC(y-1,x,z));
                    end
                end 
            end
            Diffyx(y,x,z) = diff_1x;
            Diffyy(y,x,z) = diff_1y;
            Diffyz(y,x,z) = diff_1z;
        end
    end
end

Diffyx(isnan(Diffyx)) = 0;
Diffyy(isnan(Diffyy)) = 0;
Diffyz(isnan(Diffyz)) = 0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%end of file