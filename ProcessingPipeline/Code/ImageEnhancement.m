% BreastROIHR_SUB(BreastROIHR_SUB<0) = 0;
% img = BreastROIHR_SUB;
function [I_enh] = ImageEnhancement(img0, k)
img = img0;

[Ydim, Xdim, Zdim] =  size(img);
img(img <= 0) = NaN;

Max_perslice = max(max(img));
img = 255 * img ./ repmat(Max_perslice, [Ydim, Xdim, 1]);

%img = 255 * img / max(img(:));
%% lavg
kavg_x = k;
kavg_y = k;
kavg_z = k;

lavg = zeros(Ydim, Xdim, Zdim);
lmin = zeros(Ydim, Xdim, Zdim);
lmax = zeros(Ydim, Xdim, Zdim);

for z = 1:Zdim
%     z
    for y = 1:Ydim
        for x = 1:Xdim
            Y = max(1, y-kavg_y):min(Ydim, y+kavg_y);
            X = max(1, x-kavg_x):min(Xdim, x+kavg_x);
            Z = max(1, z-kavg_z):min(Zdim, z+kavg_z);
            
            neighbor = img(Y, X, Z);
            
            lavg(y,x,z) = mean(neighbor(:), 'omitnan'); 
            lmin(y,x,z) = min(neighbor(:));
            lmax(y,x,z) = max(neighbor(:));
        end
    end
end

%% Gaussain filter
lmax(isnan(lmax)) = img0(isnan(lmax));
lmin(isnan(lmin)) = img0(isnan(lmin));
lavg(isnan(lavg)) = img0(isnan(lavg));

kernel = 2;
lmax = imgaussfilt3(lmax, kernel);
lmin = imgaussfilt3(lmin, kernel);
lavg = imgaussfilt3(lavg, kernel);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% range transfer
Range_orig = abs(lmax - lmin);
omiga0 = 55;

Range_trans = (Range_orig<=omiga0).*(omiga0 - sqrt(omiga0^2 - Range_orig.^2)) ...
             +(Range_orig>omiga0) .*(omiga0 + sqrt((255-omiga0)^2 - (255-Range_orig).^2));
         
I_new = Range_trans .* (img - lmin) ./ Range_orig;
A_new = Range_trans .* (lavg - lmin) ./ Range_orig;

%%
alpha = (A_new - I_new) / 128;
a = alpha ./ (2 * Range_trans);
b = I_new .* alpha ./ Range_trans - alpha - 1;
c = I_new.^2 .* alpha ./ (2 * Range_trans) - alpha .* I_new + I_new;
I_enh = lmin + (-b - sqrt(b.^2 -4*a.*c)) ./ (2*a);

%% save
%{
for i = 1:Zdim

% imagesc(img1)
% imshow(uint8(round(255*(img1-scale_down)/scale)))
% pause(0.05)
imwrite(uint8(BreastROIHR_SUB(:,:,i)), ['Img_slice', num2str(i), '.jpg'])

%imwrite(uint8(round(255*(Contrast(:,:,i)))), ['imadjust_slice', num2str(i), '.jpg'])
end
%}
      

