function ADC_map = CalculateADC(DWI)

DWI_Data = double(DWI.Data);
[ysize, xsize, zsize, numbval] = size(DWI_Data);
DWI.bval = [0 100 600 800]; %%% CHANGE HERE %%%
bval     = DWI.bval;

%% ADC fitting
numvox = ysize*xsize*zsize;
DWI_logData = log(reshape(DWI_Data,[numvox,numbval]));
R = zeros(numvox,1);
B = zeros(numvox,1);
ADC = zeros(numvox,1);
Fitting = zeros(size(DWI_logData));

for ii=1:numvox
    if ~mod(ii, 10000)
        ii
    end
    
    logData_ind = DWI_logData(ii,:);
    
    
    if sum(diff(logData_ind)>0)>0 || sum(isinf(logData_ind))>0
        R(ii) = 0; ADC(ii) = -1; B(ii) = 0;
        Fitting(ii,:) = NaN * ones(size(logData_ind));
    else
        [R(ii),ADC(ii),B(ii)] = regression(-bval,logData_ind);
        Fitting(ii,:) = ADC(ii)*(-bval) + B(ii);
    end
    
end

ADC_map = reshape(ADC, [ysize, xsize, zsize]);


ADC_map(ADC_map < 0) = 0;
ADC_map = ADC_map * 1e6;

%                




