matpath = '/Users/chengyuewu/Desktop/study_work/Private/MDACC_temp/Presentation/ISBI2025/Demo/ProcessingPipelineClean/Data/';
CaseNum = '104268T2';

load([matpath, CaseNum, '_DCE_scandata.mat'])
DCEsubtraction = max(DCE.Data,[],4) - DCE.Data(:,:,:,1);

%%
% for i = [50:(-1):1, 50:size(DCEsubtraction,3)]
%     imagesc(DCEsubtraction(:,:,i))
%     title(i)
%     pause
% end

    % TumorSlices  = 40:60;
    % TumorSlices  = 42:64;
    TumorSlices  = 43:57;

%% 1. Draw a approximate margin of tumor
[ysize, xsize, zsize, tsize] = size(DCE.Data);
TumorMargin = zeros(ysize, xsize, zsize);
    
    figure
    imagesc(max(DCEsubtraction(:,:,TumorSlices),[],3));
    
    temp = roipoly();
    TumorMargin(:,:,TumorSlices) = repmat(temp,[1,1,length(TumorSlices)]);


% % %%
% % for i = range_hs(round(length(range_hs)/2)-1):(-1):range_hs(1) %1:zsize
% % 
% %     subtraction = (DYNHR_br(:,:,end,i) - DYNHR_br(:,:,1,i));
% %     imagesc(subtraction)
% %     colorbar
% %     colormap jet
% %     caxis([0 max(subtraction(:))]);
% %     title(i)
% %     temp = roipoly();
% %     TumorMargin(:,:,i) = temp;
% % 
% % end

%% 2. Enhancement by dividing the post-contrast intensity by the pre-contrast
Ienh = DCE.Data .* repmat(TumorMargin,[1, 1, 1, tsize]);
Ienh = Ienh - Ienh(:,:,:,1);
Ienh = reshape(Ienh,[xsize*ysize*zsize, tsize]);

Condition = sum(Ienh,2);
Condition(isnan(Condition)) = 0;

vox = find(Condition);
X   = Ienh(vox,:);

%% 3. 2-category FCM
[V,U,~] = FCMClust(X,2);

%% 4. binarization of lesion membership
threshold = 0.3;
LMM = zeros(ysize*xsize*zsize,1);
if sum(V(1,:).^2) > sum(V(2,:).^2)
    LMM(vox) = U(1,:);
else 
    LMM(vox) = U(2,:);
end

LMM = reshape(LMM,[ysize xsize zsize]);


% % BW_lesion = (imgaussfilt3(LMM,1)>threshold);
BW_lesion = (LMM > threshold);

%% 5. reduce the suprious structure
BW_lesion_ei = bwareaopen(BW_lesion, 27);

%% 6. hole-filling operation
% BW_lesion_hf = zeros(size(BW_lesion_ei));
% for z = 1:size(BW_lesion_ei,3)
%     BW_lesion_hf(:,:,z) = imfill(BW_lesion_ei(:,:,z),'holes');
% end

%% final
% TumorROI = BW_lesion_hf;
TumorROI = BW_lesion_ei;

%% save results
save([matpath CaseNum,'_ParaMap_DCEmask_scandata.mat'], 'TumorROI','TumorMargin','U','V','threshold','LMM')

%%
close all
figure
for slice = round(median(TumorSlices)) %range_hs

    subtraction = DCEsubtraction(:,:,slice);

    subplot(1,3,1)
    [Y,X] = find(edge(TumorMargin(:,:,slice)));
    imagesc(subtraction) %TumorMargin(:,:,slice).*
    hold on
    plot(X, Y, 'r.')
    caxis([0 max(subtraction(:))]);
    colorbar
    axis off
    
    subplot(1,3,2)
    imagesc(LMM(:,:,slice))
    caxis([0 1])
    colorbar
    axis off
    
    subplot(1,3,3)
    imagesc(TumorROI(:,:,slice).*subtraction)
    caxis([0 max(subtraction(:))]);
    colorbar
    axis off
    
    % colormap(c)
    % pause
end


%%
%
testvol = TumorROI;
hold on
h_final = figure;
col=[0.8 0.8 0.8];

hiso = patch(isosurface(testvol,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(testvol,0),'FaceColor',col,'EdgeColor','none');
% axis equal;
axis on; grid on
lighting phong;
isonormals(testvol,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
view([0.5 0.8 0.5])
pause(1)

%}


