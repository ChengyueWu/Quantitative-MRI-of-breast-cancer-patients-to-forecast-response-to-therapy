function [imagesizes,scan1,scan2,scan3] = GenerateScanStructs_postReg_ISPY2(currentpath,subname)    

path = ([currentpath '/' subname '/IntervisitRegistered/RegFiles/']);

scan1 = struct;
scan2 = struct;
scan3 = struct;

%For all four scans, load different data files
%Scans 1 and 3 files named differently from Scan 2 due to registration
for v = 1:3
    fprintf('Loading Scan #%1i Images .... \n',v);
    %%% CHANGE HERE %%%
    filename = [path subname 'T' num2str(v-1) '_Reg_Manual_Weight50.mat'];
%     filename = [path subname 'T' num2str(v-1) '_Reg.mat'];
    load(filename);
    
    %Save to the scan structure
    eval(['scan' num2str(v) '.avgdce = avgdce;']);
    eval(['scan' num2str(v) '.enhanced = anatomical;']);
    eval(['scan' num2str(v) '.roi = roi;']);
    eval(['scan' num2str(v) '.adc = adc/(10^6);']); 
            
end

imagesizes = size(scan1.adc);