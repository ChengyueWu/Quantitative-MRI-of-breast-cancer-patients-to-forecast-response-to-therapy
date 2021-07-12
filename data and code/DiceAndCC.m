%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Function to calculate overlap and correlation between two matrices

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [Dice,PCC,CCC]=DiceAndCC(A,B)
% Inputs:
%      A  :     matrix of values                    double  (sy,sx*sz)
%      B  :     matrix of values                    double  (sy,sx*sz)
% Outputs:
%      Dice     :     Dice coefficient              double  (1x1)
%      PCC      :     Pearson CC                    double  (1x1)
%      CCC      :     Concordance CC                double  (1x1)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Calculate correlation between the two matrices
%Pearson's correlation coefficient
PCC = corr2(A,B);
%Concordance correlation coefficient
CCC = ccc_barnes2(reshape(A,[],1),reshape(B,[],1));

%Calculations for the Dice coefficient for overlap
A(A>1) = 1;
A(A<1) = 0;
A(isnan(A)) = 0;
B(B>1) = 1;
B(B<1) = 0;
B(isnan(B)) = 0;

NumAllElements = A + B;

intersectTC = A + B;
intersectTC(intersectTC<=1) = 0;
intersectTC(intersectTC>0) = 1;

Dice = 2*sum(intersectTC)/sum(NumAllElements);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end of file