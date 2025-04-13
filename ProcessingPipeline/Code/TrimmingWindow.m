function [X,Y] = TrimmingWindow(scan2)   


%TRIMMING DATA MATRICES TO WINDOW SIZE

%Defining a tight window that includes entire breast mask across all slices
%to trim matrices to this square region only
TotalMask = sum(scan2.maskbreast,3);
[rows,columns] = find(TotalMask>0);
xmin = min(columns);
xmax = max(columns);
ymin = min(rows);
ymax = max(rows);

%Finding longest length
xdiff = xmax - xmin;
ydiff = ymax - ymin;
n = round(max([xdiff ydiff])) + 1;
%The parameter n is the square window size for the patient


%Adjust smaller side to make domain square and generate gridspace
if ydiff > xdiff
    extra = 0;
    remainder = ydiff - xdiff;
    xmin = xmin - ceil(remainder/2);
    if xmin < 1
        extra = 1 - xmin;
        xmin = 1;
    elseif xmin == 0
        extra = 1;
        xmin = 1;
    end
    xmax = xmax + (remainder - ceil(remainder/2)) + extra;
    if xmax > size(scan2.maskbreast,2)
       extra = xmax - size(scan2.maskbreast,2);
       xmax = size(scan2.maskbreast,2);
       xmin = xmin - extra;
    end
elseif xdiff > ydiff 
    extra = 0;
    remainder = xdiff - ydiff;
    ymin = ymin - ceil(remainder/2);
    if ymin < 1
        extra = 1 - ymin;
        ymin = 1;
    elseif ymin == 0
        extra = 1;
        ymin = 1;
    end
    ymax = ymax + (remainder - ceil(remainder/2))+extra;
    if ymax > size(scan2.maskbreast,1)
       extra = ymax - size(scan2.maskbreast,1);
       ymax = size(scan2.maskbreast,1);
       ymin = ymin - extra;
    end
end    

%Creating a meshgrid for the position of the region of interest across all
%the scan data
[X,Y] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));

[d1, d2, d3] = size(scan2.maskbreast);
X(X < 1) = 1; X(X > d2) = d2;
Y(Y < 1) = 1; Y(Y > d1) = d1;



