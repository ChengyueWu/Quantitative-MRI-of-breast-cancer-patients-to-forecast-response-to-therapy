function [scan1,scan2,scan3,cellr,cellVolume,voxVolume,CarryCap,ADCw,ADCmin] = CalculateNTCs(scan1,scan2,scan3,ADCdims)   
%CALCULATE TUMOR CELL NUMBERS

%Define approximate cell radius and volume
cellr = 0.010; %millimeters, mm
cellVolume = (4/3)*pi*(cellr^3);
%For DCE data dimensions (which the ADC has been converted to)
voxVolume = ADCdims(1)*ADCdims(2)*ADCdims(3);  %mm resolution of imaging voxels
%Calculate carrying capacity with approximately packing density
CarryCap = (0.7405*voxVolume)/cellVolume;
%Mean water ADC (Nkiruka)
ADCw = 3e-3;  

%First identify slices with tumor ROIs and zeros in the ROI
if isfield(scan1,'roi')
    sl = [];
    for i = 1:size(scan1.roi,3)
        if max(max(scan1.roi(:,:,i)))>0
            sl = [sl i];
        end
    end
    scan1.tumorbearingslices = sl;

    roinan = scan1.roi;
    roinan(roinan == 0) = NaN;
    if numel(find(scan1.adc.*roinan==0)) > 0
        %Had to update this section from the UT version because we cannot
        %assume the data is square in the XY plane
        a = find(scan1.adc.*roinan==0);
        [rows,colms,slices] = ind2sub(size(scan1.adc),find(scan1.adc.*roinan==0));
        for gg = 1:numel(a)
            neighbors = zeros(1,8);
            neighbors(1) = scan1.adc(rows(gg)-1,colms(gg)-1,slices(gg));
            neighbors(2) = scan1.adc(rows(gg)+1,colms(gg)+1,slices(gg));
            neighbors(3) = scan1.adc(rows(gg),colms(gg)+1,slices(gg));
            neighbors(4) = scan1.adc(rows(gg),colms(gg)-1,slices(gg));
            neighbors(5) = scan1.adc(rows(gg)-1,colms(gg),slices(gg));
            neighbors(6) = scan1.adc(rows(gg)+1,colms(gg),slices(gg));
            neighbors(7) = scan1.adc(rows(gg)-1,colms(gg)+1,slices(gg));
            neighbors(8) = scan1.adc(rows(gg)+1,colms(gg)-1,slices(gg));
            scan1.adc(rows(gg),colms(gg),slices(gg)) = mean(neighbors);
        end
    end
end

if isfield(scan2,'roi')
    sl = [];
    for i = 1:size(scan2.roi,3)
        if max(max(scan2.roi(:,:,i)))>0
            sl = [sl i];
        end
    end
    scan2.tumorbearingslices = sl;

    roinan = scan2.roi;
    roinan(roinan == 0) = NaN;
    if numel(find(scan2.adc.*roinan==0)) > 0
        a = find(scan2.adc.*roinan==0);
        [rows,colms,slices] = ind2sub(size(scan2.adc),find(scan2.adc.*roinan==0));
        for gg = 1:numel(a)
            neighbors = zeros(1,8);
            neighbors(1) = scan2.adc(rows(gg)-1,colms(gg)-1,slices(gg));
            neighbors(2) = scan2.adc(rows(gg)+1,colms(gg)+1,slices(gg));
            neighbors(3) = scan2.adc(rows(gg),colms(gg)+1,slices(gg));
            neighbors(4) = scan2.adc(rows(gg),colms(gg)-1,slices(gg));
            neighbors(5) = scan2.adc(rows(gg)-1,colms(gg),slices(gg));
            neighbors(6) = scan2.adc(rows(gg)+1,colms(gg),slices(gg));
            neighbors(7) = scan2.adc(rows(gg)-1,colms(gg)+1,slices(gg));
            neighbors(8) = scan2.adc(rows(gg)+1,colms(gg)-1,slices(gg));
            scan2.adc(rows(gg),colms(gg),slices(gg)) = mean(neighbors);
        end
    end
end

if isfield(scan3,'roi')
    sl = [];
    for i = 1:size(scan3.roi,3)
        if max(max(scan3.roi(:,:,i)))>0
            sl = [sl i];
        end
    end
    scan3.tumorbearingslices = sl;

    roinan = scan3.roi;
    roinan(roinan == 0) = NaN;
    if numel(find(scan3.adc.*roinan==0)) > 0
        a = find(scan3.adc.*roinan==0);
        [rows,colms,slices] = ind2sub(size(scan3.adc),find(scan3.adc.*roinan==0));
        for gg = 1:numel(a)
            neighbors = zeros(1,8);
            neighbors(1) = scan3.adc(rows(gg)-1,colms(gg)-1,slices(gg));
            neighbors(2) = scan3.adc(rows(gg)+1,colms(gg)+1,slices(gg));
            neighbors(3) = scan3.adc(rows(gg),colms(gg)+1,slices(gg));
            neighbors(4) = scan3.adc(rows(gg),colms(gg)-1,slices(gg));
            neighbors(5) = scan3.adc(rows(gg)-1,colms(gg),slices(gg));
            neighbors(6) = scan3.adc(rows(gg)+1,colms(gg),slices(gg));
            neighbors(7) = scan3.adc(rows(gg)-1,colms(gg)+1,slices(gg));
            neighbors(8) = scan3.adc(rows(gg)+1,colms(gg)-1,slices(gg));
            scan3.adc(rows(gg),colms(gg),slices(gg)) = mean(neighbors);
        end
    end
end    

ADC1 = scan1.adc;  
ADC2 = scan2.adc;  
ADC3 = scan3.adc;

%Eliminating any voxels that had a negative ADC value (may still be
%negative due to filter) and multiply by tumor roi.
ADC1f = ADC1.*scan1.roi;
ADC2f = ADC2.*scan2.roi; 
ADC3f = ADC3.*scan3.roi;

%Calculate mean and std of ADC values > 0
ADC1mv = mean(ADC1f(find(ADC1f>0))); 
ADC2mv = mean(ADC2f(find(ADC2f>0))); 
std1 = std(ADC1f(find(ADC1f>0))); 
std2 = std(ADC2f(find(ADC2f>0))); 
%Must check for later scans if there is tumor
if max(max(max(scan3.roi))) == 0 
    ADC3mv = 1; 
    std3 = 0; 
else
    ADC3mv = mean(ADC3f(find(ADC3f>0))); 
    std3 = std(ADC3f(find(ADC3f>0)));
end

%Calculating an absolute minimum and using 95% confidence: mean - 1.96*std
min1a = min(min(min(ADC1f(find(ADC1f>0)))));
min1c = ADC1mv-(1.96*std1);
if min1c<0 
    min1c = 1; 
end
min2a = min(min(min(ADC2f(find(ADC2f>0)))));
min2c = ADC2mv-(1.96*std2);
if min2c<0
    min2c = 1; 
end
min3a = min(min(min(ADC3f(find(ADC3f>0)))));
min3c = ADC3mv-(1.96*std3);
if min3c<0
    min3c = 1;
end

%Find minimum over all days to calculate NTC
%Absolute minimum
ADCmina = min([min1a min2a min3a]);
%Calculate minimum using 95% confidence
ADCminc = min([min1c min2c min3c]);
ADCmin  = min([ADCminc ADCmina]);

%Using ADC minimum, water ADC, and approximate carrying capacity, calculate
%the number of tumor cells per voxel
N1 = CarryCap*(ADCw-ADC1)/(ADCw-ADCmin); 
N2 = CarryCap*(ADCw-ADC2)/(ADCw-ADCmin); 
N3 = CarryCap*(ADCw-ADC3)/(ADCw-ADCmin); 

%Only save positive values within the tumor ROI !!! And with valid ADC value
scan1.NTC = N1.*(N1>0).*scan1.roi .* (ADC1>0); 
scan2.NTC = N2.*(N2>0).*scan2.roi .* (ADC2>0);
scan3.NTC = N3.*(N3>0).*scan3.roi .* (ADC3>0);

    