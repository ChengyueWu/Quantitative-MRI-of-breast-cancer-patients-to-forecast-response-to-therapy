function [Xlist, Ylist, Zlist, EdgeXY, Xmap, Ymap] = LocationExtract_ISPY2(InputSet,flipped) 

    y0   = double(InputSet.pos(1,2,1));
    x0   = double(InputSet.pos(1,1,1));
    res  = double(InputSet.resolution(1));
    FOV_x= double(InputSet.columns * res);
    FOV_y= double(InputSet.row * res);

    if flipped == 1
        Zlist = double(InputSet.pos(:,3,1));
        Xlist = x0 - (0:double(InputSet.columns-1))     * res;
        Ylist = y0 - (0:double(InputSet.row-1)) * res;
        [Ymap, Xmap] = meshgrid(Ylist, Xlist);
        
        EdgeXY = [x0+res/2        , y0+res/2; ...
                  x0+res/2 - FOV_x, y0+res/2; ...
                  x0+res/2 - FOV_x, y0+res/2 - FOV_y; ...
                  x0+res/2        , y0+res/2 - FOV_y; ...
                  ];
    else
        Zlist = double(InputSet.pos(:,3,1));
        Xlist = x0 + (0:double(InputSet.columns-1))     * res;
        Ylist = y0 + (0:double(InputSet.row-1)) * res;
        [Ymap, Xmap] = meshgrid(Ylist, Xlist);
        
        EdgeXY = [x0+res/2        , y0+res/2; ...
                  x0+res/2 + FOV_x, y0+res/2; ...
                  x0+res/2 + FOV_x, y0+res/2 + FOV_y; ...
                  x0+res/2        , y0+res/2 + FOV_y; ...
                  ];
    end