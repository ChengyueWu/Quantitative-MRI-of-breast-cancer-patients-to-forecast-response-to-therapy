%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Function for identifying the boundaries of the breast ROI

%Recieves a 0 and 1 matrix indicating the breast region versus not
%Searches counter clockwise to identify the boundary as well as the
%boundary slices.
%Defines the boundary types using the formatting required in the 
%forwardsolve function.

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [BCF] = Boundaries3DUT(n,Zbw)
% Inputs:
%      Zbw     :     Breast mask                     logical (sy,sx,sz)
%      n       :     grid size                       double  (1x1)
% Outputs:
%      BCF     :     Boundary matrix                 double  (sy,sx,sz)
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

%Define boundary matrix the same size as breast mask
BCF = zeros(size(Zbw));
gridsize = [size(Zbw,1),size(Zbw,2)];

%Zeroing out edges 
Zbw(end,:,:) = 0;
Zbw(1,:,:) = 0;
Zbw(:,end,:) = 0;
Zbw(:,1,:) = 0;

%loop through each slice to define the boundary
for jj = 1:size(Zbw,3)
    clearvars -EXCEPT Zbw BCF n gridsize jj tvec timesteps X Y dx dy
    %``````````````````````````````````````````````````````````````````````````
    %Holding matrix for the breast mask slice
    ZBW = Zbw(:,:,jj);
    ZBW(isnan(ZBW)) = 0;
    
    %Default boundary values are 5 (inner voxel) or zero (outside domain)
    BCF(:,:,jj) = 5*ZBW;

    %Make a boundary tracing
    tracing = bwperim(ZBW);
    %Convert from logical to numerical
    tracing = double(tracing);

    %For the edges remove the value 5
    BCF(:,:,jj) = BCF(:,:,jj) - 5*tracing;

    %Redefining the shaoe Zbw
    ZBW = reshape(ZBW',[],1);
    temp = zeros(size(ZBW));

    %Flags to catch ends of the breast mask
    flag = 0;
    newflag = 0;
    count = 1;
    newcount = 1;
    while(newflag==0)
        flag = 0;
        count = 1;
        while(flag==0)
            test = tracing(count,newcount);
            if(test==1)
                start = (count-1)*gridsize(1) + newcount;
                flag = 1;
                newflag = 1;
            elseif(count==n)
                flag = 1;
            end    
            count = count + 1; 
        end    
        newcount = newcount + 1;
    end

    %Find each of the voxels in the boundary tracing
    bnodes = zeros(1,size(find(tracing==1),1)+2);
    bnodesx = zeros(1,size(find(tracing==1),1)+2);
    bnodesy = zeros(1,size(find(tracing==1),1)+2);
    bnodes(3) = start;
    bnodesx(3) = newcount - 1;
    bnodesy(3) = count - 1;

    %Counter-clockwise algorithm to determine what type of edge the voxel
    %in the tracing is
    for ii = 3:size(find(tracing==1),1)+1
        rowcur = floor(bnodes(ii)/gridsize(2))+1;
        colcur = mod(bnodes(ii),gridsize(1));
        for position = 1:8
            if(position==1&&(rowcur-1<1||colcur+1>n))
                 continue;   
            end
            if(position==2&&colcur+1>n)
                 continue;    
            end         
            if(position==3&&(rowcur+1>n||colcur+1>n))
                 continue;   
            end
            if(position==4&&rowcur+1>n)
                 continue;    
            end
            if(position==5&&(rowcur+1>n||colcur-1<1))
                 continue;   
            end        
            if(position==6&&colcur-1<1)
                 continue;     
            end
            if(position==7&&(rowcur-1<1||colcur-1<1))
                 continue;    
            end          
            if(position==8&&rowcur-1<1)
                 continue;    
            end
            if((position==1)&&(tracing(rowcur-1,colcur+1)~=0)&&((((rowcur-1)-1)*gridsize(2)+colcur+1)~=bnodes(ii-1))&&((((rowcur-1)-1)*gridsize(2)+colcur+1)~=bnodes(ii-2)))
            %down 1, right 1
                bnodes(ii+1) = ((rowcur-1)-1)*gridsize(2) + colcur + 1;
                bnodesx(ii+1) = colcur+1;
                bnodesy(ii+1) = rowcur-1;
                break
            elseif((position==2)&&(tracing(rowcur,colcur+1)~=0)&&((((rowcur)-1)*gridsize(2)+colcur+1)~=bnodes(ii-1))&&((((rowcur)-1)*gridsize(2)+colcur+1)~=bnodes(ii-2)))
            %down 0, right 1
                bnodes(ii+1) = ((rowcur)-1)*gridsize(2) + colcur + 1;    
                bnodesx(ii+1) = colcur+1;
                bnodesy(ii+1) = rowcur;  
                break
            elseif((position==3)&&(tracing(rowcur+1,colcur+1)~=0)&&((((rowcur+1)-1)*gridsize(2)+colcur+1)~=bnodes(ii-1))&&((((rowcur+1)-1)*gridsize(2)+colcur+1)~=bnodes(ii-2)))
            %up 1, right 1
                bnodes(ii+1) = ((rowcur+1)-1)*gridsize(2) + colcur + 1;
                bnodesx(ii+1) = colcur+1;
                bnodesy(ii+1) = rowcur+1;
                break
            elseif((position==4)&&(tracing(rowcur+1,colcur)~=0)&&((((rowcur+1)-1)*gridsize(2)+colcur)~=bnodes(ii-1))&&((((rowcur+1)-1)*gridsize(2)+colcur)~=bnodes(ii-2)))
            %up 1, right 0    
                bnodes(ii+1) = ((rowcur+1)-1)*gridsize(2) + colcur;
                bnodesx(ii+1) = colcur;
                bnodesy(ii+1) = rowcur+1;
                break
            elseif((position==5)&&(tracing(rowcur+1,colcur-1)~=0)&&((((rowcur+1)-1)*gridsize(2)+colcur-1)~=bnodes(ii-1))&&((((rowcur+1)-1)*gridsize(2)+colcur-1)~=bnodes(ii-2)))
            %up 1, left 1    
                bnodes(ii+1) = ((rowcur+1)-1)*gridsize(2) + colcur - 1;
                bnodesx(ii+1) = colcur-1;
                bnodesy(ii+1) = rowcur+1;
                break
            elseif((position==6)&&(tracing(rowcur,colcur-1)~=0)&&((((rowcur)-1)*gridsize(2)+colcur-1)~=bnodes(ii-1))&&((((rowcur)-1)*gridsize(2)+colcur-1)~=bnodes(ii-2)))
            %up 0, left 1  
                bnodes(ii+1) = ((rowcur)-1)*gridsize(2) + colcur - 1;
                bnodesx(ii+1) = colcur-1;
                bnodesy(ii+1) = rowcur;
                break
            elseif((position==7)&&(tracing(rowcur-1,colcur-1)~=0)&&((((rowcur-1)-1)*gridsize(2)+colcur-1)~=bnodes(ii-1))&&((((rowcur-1)-1)*gridsize(2)+colcur-1)~=bnodes(ii-2)))
            %down 1, left 1    
                bnodes(ii+1) = ((rowcur-1)-1)*gridsize(2) + colcur - 1; 
                bnodesx(ii+1) = colcur-1;
                bnodesy(ii+1) = rowcur-1;
                break
            elseif((position==8)&&(tracing(rowcur-1,colcur)~=0)&&((((rowcur-1)-1)*gridsize(2)+colcur)~=bnodes(ii-1))&&((((rowcur-1)-1)*gridsize(2)+colcur)~=bnodes(ii-2)))
            %down 1, left 0
                bnodes(ii+1) = ((rowcur-1)-1)*gridsize(2) + colcur;
                bnodesx(ii+1) = colcur;
                bnodesy(ii+1) = rowcur-1;
                break
            end    
        end   
    end
    %removing placeholders
    bnodes(1:2) = [];
    bnodes(end) = [];
    bnodesx(1:2) = [];
    bnodesx(end) = [];
    bnodesy(1:2) = [];
    bnodesy(end) = [];

    %Identifying portions of the boundary for the boundary condition
    %implementation
    tops = zeros(1,2);
    bums = zeros(1,2);
    lefs = zeros(1,2);
    righs = zeros(1,2);
    blcorns = zeros(1,2);
    brcorns = zeros(1,2);
    trcorns = zeros(1,2);
    tlcorns = zeros(1,2);
    for bb = 1:size(bnodes,2)
        if(bb==1)
            before = bnodes(end);
            after = bnodes(bb+1);
        elseif(bb==size(bnodes,2))
            before = bnodes(bb-1);
            after = bnodes(1);
        else
            before = bnodes(bb-1);
            after = bnodes(bb+1);
        end

        %Corners-Bottom Left
        if(((before==(bnodes(bb)-1))||(before==(bnodes(bb)-gridsize(1)-1)))...
             &&((after==(bnodes(bb)+gridsize(1)+1))||(after==(bnodes(bb)+gridsize(1)))))            
            brcorns(end+1,1) = bnodes(bb);
            brcorns(end,2) = bb;    
        %Corners-Bottom Right
        elseif(((before==(bnodes(bb)-gridsize(1)))||(before==(bnodes(bb)-gridsize(1)+1)))...
             &&((after==(bnodes(bb)-1))||(after==(bnodes(bb)+gridsize(1)-1))))
            trcorns(end+1,1) = bnodes(bb);
            trcorns(end,2) = bb;    
        %Corners-Top Right
        elseif(((before==(bnodes(bb)+1))||(before==(bnodes(bb)+gridsize(1)+1)))...
             &&((after==(bnodes(bb)-gridsize(1)))||(after==(bnodes(bb)-gridsize(1)-1))))
            tlcorns(end+1,1) = bnodes(bb);
            tlcorns(end,2) = bb;    
        %Corners-Top Left    
        elseif(((before==(bnodes(bb)+gridsize(1)))||(before==(bnodes(bb)+gridsize(1)-1)))...
            &&((after==(bnodes(bb)+1))||(after==(bnodes(bb)-gridsize(1)+1))))       
            blcorns(end+1,1) = bnodes(bb);
            blcorns(end,2) = bb;
        %Horizontal Edges-top
        elseif((after==(bnodes(bb)-gridsize(1)))||(before==(bnodes(bb)+gridsize(1))))    
            lefs(end+1,1) = bnodes(bb);
            lefs(end,2) = bb;
        %Horizontal Edges-Bottom
        elseif((before==(bnodes(bb)-gridsize(1)))||(after==(bnodes(bb)+gridsize(1))))        
            righs(end+1,1) = bnodes(bb);
            righs(end,2) = bb;
        %Vertical Edges-Left    
        elseif((before==(bnodes(bb)-1))||(after==(bnodes(bb)+1)))        
            bums(end+1,1) = bnodes(bb);
            bums(end,2) = bb;
        %Vertical Edges-Right    
        elseif(after==(bnodes(bb)-1))||(before==(bnodes(bb)+1))  
            tops(end+1,1) = bnodes(bb);
            tops(end,2) = bb;
        end  
    end    
    tops(1,:) = [];
    bums(1,:) = [];
    lefs(1,:) = [];
    righs(1,:) = [];
    blcorns(1,:) = [];
    brcorns(1,:) = [];
    trcorns(1,:) = [];
    tlcorns(1,:) = [];

    downers = sortrows([blcorns;bums;brcorns],2);
    rights = sortrows([brcorns;righs;trcorns],2);
    uppers = sortrows([trcorns;tops;tlcorns],2);
    lefts = sortrows([tlcorns;lefs;blcorns],2);

    temp(downers(:,1)) = temp(downers(:,1)) + 1;
    temp(uppers(:,1)) = temp(uppers(:,1)) + 2;
    temp(lefts(:,1)) = temp(lefts(:,1)) + 20;
    temp(rights(:,1)) = temp(rights(:,1)) + 10;

    %update BCF
    BCF(:,:,jj) = BCF(:,:,jj) + reshape(temp,gridsize)';
end

BCF(1,:,:) = BCF(1,:,:) + 1;
BCF(end,:,:) = BCF(end,:,:) + 2;
BCF(:,1,:) = BCF(:,1,:) + 20;
BCF(:,end,:) = BCF(:,end,:) + 10;

%Defining what slices might be above or below a given slice to determine
%z direction boundaries
for ii = 1:size(Zbw,3)
    ZBW = Zbw(:,:,ii);
    ZBW(isnan(ZBW)) = 0;
    temp = zeros(size(ZBW));
    
    if ii == 1
        above = Zbw(:,:,ii+1);
        above(isnan(above)) = 0;
        topper = ZBW - above;
        a = find(topper == 1);
        b = find(ZBW==1);
    elseif ii == size(Zbw,3)
        below = Zbw(:,:,ii-1);
        below(isnan(below)) = 0;
        bummer = ZBW - below;
        a = find(ZBW==1);
        b = find(bummer == 1);
    else
        above = Zbw(:,:,ii+1);
        above(isnan(above)) = 0;
        topper = ZBW - above;
        a = find(topper == 1); 
        below = Zbw(:,:,ii-1);
        below(isnan(below)) = 0;
        bummer = ZBW - below;
        b = find(bummer == 1);
    end        
    
    temp(a) = temp(a) + 100;
    temp(b) = temp(b) + 200;
    
    BCF(:,:,ii) = BCF(:,:,ii) + temp;
    
end

BCF(end,:,:) = 0;
BCF(1,:,:) = 0;
BCF(:,end,:) = 0;
BCF(:,1,:) = 0;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end of file