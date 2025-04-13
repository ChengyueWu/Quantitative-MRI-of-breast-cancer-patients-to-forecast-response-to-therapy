function [scan1,scan2,scan3] = SmoothingOutMaps(scan1,scan2,scan3)

A = fieldnames(scan1);
B = fieldnames(scan2);
C = fieldnames(scan3);

a = find(strcmp(A,'roi'));
b = find(strcmp(A,'maskbreast'));
c = find(strcmp(A,'maskfibro'));
d = find(strcmp(A,'maskadipo'));
% e = find(strcmp(A,'enhanced'));
% f = find(strcmp(A,'avgdce'));
g = find(strcmp(A,'NTC'));
h = find(strcmp(A,'tumorbearingslices'));

which = 1:size(A,1);
which([a,b,c,d,g,h]) = []; %e,f,

A = A(which);
B = B(which);
C = C(which);

for jj = 1:size(A,1)
    % Smooth data maps
    hh = fspecial('gaussian',[3 3],1.5);
    
    if length(size(eval(['scan1.' A{jj}]))) == 3
    %     eval(['temp1 = zeros(size(scan1.', A{jj}, '));')
        for ii = 1:size(scan1.roi,3)
            eval(['temp1(:,:,ii) = imfilter(scan1.' A{jj} '(:,:,ii),hh);'])
            eval(['temp2(:,:,ii) = imfilter(scan2.' B{jj} '(:,:,ii),hh);'])
            eval(['temp3(:,:,ii) = imfilter(scan3.' C{jj} '(:,:,ii),hh);'])
        end

        temp1(temp1<0) = 0;
        temp2(temp2<0) = 0;
        temp3(temp3<0) = 0;

        if strcmp(A{jj},'adc') == 1
            temp1(temp1>0.003) = 0.003;
            temp2(temp2>0.003) = 0.003;
            temp3(temp3>0.003) = 0.003;
        end

        eval(['scan1.' A{jj} '= temp1;'])
        eval(['scan2.' B{jj} '= temp2;'])
        eval(['scan3.' C{jj} '= temp3;'])
        clear temp1 temp2 temp3
    end
end


