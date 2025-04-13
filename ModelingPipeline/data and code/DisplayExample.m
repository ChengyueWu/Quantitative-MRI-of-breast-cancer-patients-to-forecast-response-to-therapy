%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.


close all

numslices = size(TC,3);
slice = median(1:numslices);

BreastMask(BreastMask==0) = NaN;
surf(flipud(TC(:,:,slice).*BreastMask(:,:,slice)),'EdgeColor','none'); 
view(2)
xticks ''
yticks ''
grid 'off'
colormap jet
a = colorbar;
ylabel(a,'Cells','FontSize',16,'Rotation',270);
a.Label.Position(1) = 4;
a.FontName = 'Times New Roman';
caxis([0 16*10^5])
title(['Prediction for ' patientID ', Slice ' num2str(slice) ' of ' num2str(numslices) ])
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',20)

clc