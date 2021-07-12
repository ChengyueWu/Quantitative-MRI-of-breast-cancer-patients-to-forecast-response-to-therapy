%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Function to calculate the concordance correlation coefficient between two
%vectors

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [ccc] = CalcCCC(x,y)
% Inputs:
%      x    :     first vector                          double  (sy,sx*sz)
%      y    :     second vector                         double  (sy,sx*sz)
% Outputs:
%      ccc  :     concordance correlation coefficient   double  (1x1)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%number of elements in x and y vectors
numx = max(size(x));
numy = max(size(y));

%convert to column vectors if inputs are row vectors
if size(x,2)>size(x,1)
    x=x';
end
if size(y,2)>size(y,1)
    y=y';
end

%if arrays are not the same length, exit function
if numx~=numy
    fprintf('\n\nError: vectors are not the same length.\n');
    ccc = 'NAN';
    return
end
if min(size(x))>1 || min(size(y))>1
    fprintf('\n\nError: inputs are matrices. Function requires vectors.\n');
    ccc = 'NAN';
    return
end

%calculate means
xmean = sum(x)/(numx);
ymean = sum(y)/(numy);

%calculate variances and covariance
sxx = (sum((x-xmean).^2))/(numx);
syy = (sum((y-ymean).^2))/(numy);
sxy = (sum((x-xmean).*(y-ymean)))/(numx);

%calculate squared difference of means
btv = (xmean-ymean)^2;

%calculate concordance correlation coefficient
ccc = (2*sxy)/(sxx+syy+btv);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end of file
