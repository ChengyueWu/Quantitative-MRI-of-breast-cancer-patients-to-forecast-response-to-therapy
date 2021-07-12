%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Function for optimization of parameter using Levenberg-Marquardt method

%%% Authors:     Angela M. Jarrett, Chengyue Wu, Thomas E. Yankeelov
%%% Last edit:   July 12, 2021
%%% Affiliation: UT Austin
%%% Reference:   Jarrett et al., "Quantitative magnetic resonance imaging
%%%              and tumor forecasting of breast cancer patients in the community
%%%              setting", Nature Protocol.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% function [] = NtcFD3DmechOPUT1to2(pt)
% if pt<10
% patientID = ['TXO_00' num2str(pt)];
% else
% patientID = ['TXO_0' num2str(pt)];    
% end
% disp(['Patient is ', patientID])

patientID = 'testpatient';

%Designate start scan
startscan = 1;
%Designate comparison scan
endscan = 2;

%Read in patient information for image dimensions and scan times
fid = fopen([patientID '.txt']);
tline = fgetl(fid);
counter = 1;
while ischar(tline)
    if(counter == 1)
        imagedims = tline;
    elseif(counter == 2)
        times = tline;
    elseif(counter == 3)
        schedule = tline;
    elseif(counter == 4)
        slicebounds = tline;  
    elseif(counter == 5)
        GridSize = tline;    
    end    
    counter = counter + 1;    
    tline = fgetl(fid);
end
fclose(fid);
imagedims = str2num(imagedims);
times = str2num(times);
%Final time in days from intial scan to comparison scan
schedule = strsplit(schedule);
schedule = char(schedule);
scans = find(schedule=='S');
tf = sum(times(scans(startscan)+1:scans(endscan))); 
slicebounds = str2num(slicebounds);
%Specific slice top and bottom for model simulation efficiency
slicestart = slicebounds(1);
sliceend = slicebounds(2);
whichslices = slicestart:sliceend;
pslices = numel(whichslices);
n = str2num(GridSize);

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Loading data, grid, and identifying boundaries

%Load native grid
load(['NativeX_' patientID '.mat']);
load(['NativeY_' patientID '.mat']);
load(['BreastMask_' patientID '.mat']);
load(['NTC1_' patientID '.mat']);
load(['Tissues1_' patientID '.mat']);
load(['NTC2_' patientID '.mat']); 
z = NTC1(:,:,whichslices);
tissues = Tissues1(:,:,whichslices);
zf = NTC2(:,:,whichslices);
bw = BreastMask(:,:,whichslices);

%Coarsening grid for computational efficiency
dx = (max(max(X))-min(min(X)))/(str2num(GridSize)/2);
dy = (max(max(Y))-min(min(Y)))/(str2num(GridSize)/2);
[xq,yq] = meshgrid(min(min(X)):dx:max(max(X)),min(min(Y)):dy:max(max(Y)));
gridsize = size(xq);
n = gridsize(1);

%Interpolate points using cubic interpolator for each slice
for kk = 1:size(bw,3)
    Z(:,:,kk) = griddata(X,Y,z(:,:,kk),xq,yq,'cubic');
    Zf(:,:,kk) = griddata(X,Y,zf(:,:,kk),xq,yq,'cubic');
    Tissues(:,:,kk) = griddata(X,Y,tissues(:,:,kk),xq,yq,'cubic');
    BW(:,:,kk) = griddata(X,Y,bw(:,:,kk),xq,yq,'cubic');
end
%Potential interpolation corrections
Z(Z<0) = 0;
Zf(Zf<0) = 0;
Tissues(Tissues<0) = 0;
BW(BW<0) = 0;

BW = round(BW);
Tissues = round(Tissues);

clear Z1 Z2 Tissues1 Tissues2 Zbw

disp('Patient data loaded')

%Identifying boundaries
[BCF] = Boundaries3DUT(n,BW);
dx = imagedims(1);
dy = imagedims(2);
dz = imagedims(3);

disp('Boundaries identified')

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Time step is a quarter day
dt = 0.25;

%Make the time vector for the simulations
tvec = 0:dt:tf;
timesteps = size(tvec,2);

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Define carrying capacity and calculate first guess
TCvolume = 4189e-9;
packingdensity = 0.7405;
voxelvolume = prod(imagedims);
tHeta = packingdensity*voxelvolume/TCvolume;

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Defining values for mechanical coupling of tissue properties to diffusion
%Shear modulus and Hooke's law matrix
nU = 0.45; 

%Convert elements for the tumor tissue
E = 2e3*ones(size(Tissues));            %Adipose
E(Tissues==2) = 2*2e3;                  %Fibroglandular
ScaledZ = Z/tHeta;
E(ScaledZ>0.1) = 10*2e3;                %Tumor
clear ScaledZ 

%Calculate shear modulus 
G = E/2*(1-nU);
clear E

%Derivatives of G
[DiffyxG,DiffyyG,DiffyzG] = Diffy3D(dx,dy,dz,BCF,G);
DiffyG = zeros(n,n,pslices,3);
DiffyG(:,:,:,1) = DiffyxG;
DiffyG(:,:,:,2) = DiffyyG;
DiffyG(:,:,:,3) = DiffyzG;
clear DiffyxG DiffyyG DiffyzG

%Define a larger ROI to focus the proliferation map calibration
ROI = zeros(size(BW));
BigROI = zeros(size(BW));
for gg = 1:pslices
    BigROI(:,:,gg) = Z(:,:,gg) + Zf(:,:,gg);
    BigROI(BigROI>0) = 1;
    BigROI(isnan(BigROI)) = 0;
    BigROIswell = imdilate(BigROI(:,:,gg),strel('disk',3));
    BigROIswell = bwconvhull(BigROIswell);
    ROI(:,:,gg) = BigROIswell.*BW(:,:,gg);
end
which = find(ROI==1);

CarryCaps = ones(n,n,pslices)*tHeta;
CarryCaps(CarryCaps<1) = 1;

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Reshaping final time info (Zstart for the calibration)
Zf = reshape(Zf,[],1);
Zf(isnan(Zf)) = 0;

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Define parameters and best guess
%Here we sample throughout the parameter space to make a "good" first guess
%for the calibration method.

disp('First guess testing')
%Using 300 intial guesses because easily divisible by 12 and between
%200-400 
runs = 300;
%Number of processors available for parallel computing
comps = 12;

%Defining different upper and lower bounds for parameters 
uD = 10^-3;
lk = -0.5;
uk = 0.1;

%params = [DTC, growth]
lower   = [10^-6, lk];
upper   = [uD, uk];

%Using uniform random distributions
RandNums = zeros(runs,size(lower,2));
for nn = 1:runs
%   r = a + (b-a).*rand(n,1);
    RandNums(nn,:) = lower + (upper - lower).*rand(1,size(lower,2));
end

for bb = 1:runs/comps
    parfor yy = 1 + comps*(bb-1):comps + comps*(bb-1)
        %Assigning parameters based on the loop
        params = RandNums(yy,:);
        DMatrix = params(1)*ones(n,n,pslices);
        kMatrix = -10*ones(n,n,pslices);
        kMatrix(which) = params(2);
        %Simulating the model
        BigTC = NtcFDmech3DUTLogisticOnly(Z,imagedims,DMatrix,dt,...
                              kMatrix,CarryCaps,BCF,timesteps,G,DiffyG,nU);    
        %Saving results
        TCtest(:,yy) = reshape(BigTC,[],1);
        errors(yy) =  sum((Zf-TCtest(:,yy)).^2);
    end
end

%Saving the best of all the intial guesses
bestguess = find(errors==min(errors));
kMatrix = -10*ones(n,n,pslices);
kMatrix(which) = RandNums(bestguess(1),2);
params = [RandNums(bestguess(1),1),reshape(kMatrix,1,[])]';
TCfinal = TCtest(:,bestguess(1));

disp('First guess successful run')

%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%Optimization, LM method

%Initialize optimization parameters
lam = 1e-20;
maxiter = 200;
iter = 1;

%Creating holding matrices and first set of test values and adjusting
%vectors to only calibrate the points of interest in the big ROI
extras = pslices*n^2 - size(which,1);
globalparams = size(params,1) - pslices*n^2;
whichones = [1,globalparams+which'];
numPs = size(params,1) - extras;
TCtest = zeros(prod(gridsize)*pslices,numPs);
J = zeros(prod(gridsize)*pslices,numPs);
testvals = params + 0.05*(1+params);

disp('Beginning for loop for parpool batches')
for bb = 1:numPs/comps
parfor yy = 1 + comps*(bb-1):comps + comps*(bb-1)
    %Assigning parameters based on the loop
    test = params;
    test(whichones(yy)) = testvals(whichones(yy));
    diff = test(whichones(yy)) - params(whichones(yy));
    DMatrix = test(1)*ones(n,n,pslices);
    kMatrix = reshape(test(2:end),n,n,pslices);
    %Simulate the model
    BigTCtest = NtcFDmech3DUTLogisticOnly(Z,imagedims,DMatrix,dt,...
                              kMatrix,CarryCaps,BCF,timesteps,G,DiffyG,nU);    
    %Save results and generate Jacobian
    TCtest(:,yy) = reshape(BigTCtest,[],1);
    J(:,yy) = (TCtest(:,yy)-TCfinal)/diff;
    
end
end

%With the parallel loops, dividing by the number of cores does not always 
%result in a clean division of parameters
disp('Potential remaining of parameters running now')
if mod(numPs,comps)~=0
parfor yy = 1+numPs-mod(numPs,comps):numPs 
    test = params;
    test(whichones(yy)) = testvals(whichones(yy));
    diff = test(whichones(yy)) - params(whichones(yy));
    DMatrix = test(1)*ones(n,n,pslices);
    kMatrix = reshape(test(2:end),n,n,pslices);
    BigTCtest = NtcFDmech3DUTLogisticOnly(Z,imagedims,DMatrix,dt,...
                              kMatrix,CarryCaps,BCF,timesteps,G,DiffyG,nU);    
    TCtest(:,yy) = reshape(BigTCtest,[],1);
    J(:,yy) = (TCtest(:,yy)-TCfinal)/diff;
end
end

%Calculate the sum of the squared error
SSE =  sum((Zf-TCfinal).^2);

disp('Calculating Jacobian')
%Calculate the Hessian
Jt = J';
Js = (Jt)*J;
alphaval = lam*trace(Js)*SSE^2;
celldiffs = Zf-TCfinal;
H = Js + sqrt(alphaval)*eye(size(Js,1),size(Js,2));

%Rescaling the Hessian for paramters with large magnitude differences
g = Jt*celldiffs;
clear Hsc gsc deltasc
gsc = zeros(size(H,1),1); 
Hsc = zeros(size(H,1),size(H,2));
for i1 = 1:size(H,1)
    gsc(i1,1) = g(i1,1)/sqrt(H(i1,i1)); 
    for i2 = 1:size(H,2)
        Hsc(i1,i2) = H(i1,i2)/(sqrt(H(i1,i1))*sqrt(H(i2,i2))); 
    end
end

disp('Calculating next del')
%Calculation of next del step in parameter space
[L, U] = lu(Hsc);
tempsc = L\gsc;
delsc = U\tempsc;
Del = zeros(size(delsc));
for i1 = 1:size(delsc,1)
    Del(i1,1) = delsc(i1,1)/sqrt(H(i1,i1)); 
end
temp = zeros(size(params));
temp(whichones) = Del;
Del = temp;
%Check diffusion upper bound
if(params(1)+Del(1)>uD)
    Del(1) = (uD-params(1))/2;
end
%Check diffusion is positive
if(params(1)+Del(1)<0)
    Del(1) = -(params(1))/2;
end
%Do not check proliferation bounds in this code

disp('Defining new test parameters')
%Define new test parameters and associated error
newtest = params;
newtest(whichones) = params(whichones) + Del(whichones);
err = norm(Del)/norm(params(whichones));

%Calculate new test results and sum of the squared errors
DMatrix = newtest(1)*ones(n,n,pslices);
kMatrix = reshape(newtest(2:end),n,n,pslices);

BigTCnewtest = NtcFDmech3DUTLogisticOnly(Z,imagedims,DMatrix,dt,...
                              kMatrix,CarryCaps,BCF,timesteps,G,DiffyG,nU);    
TCnewtest = reshape(BigTCnewtest,[],1); 
SSEnew = sum((Zf-TCnewtest).^2);

%Calculating concordance correlation coefficient of the results of the
%newest parameter set
cccTC = ccc_barnes2(Zf,TCnewtest);

%Setting intermediate value holders for the optimization process
intermedparams = newtest;
TCintermed = TCnewtest;
SSE = SSEnew; 
paramsnew = params;
lambdaup = 0;
upcnt = 0;
testcnt = 0;
recalc = 0;

disp('Beginning optimization loop')
%Optimization loop
%Here are the determinants for the calibration process: SSE, iterations,
%and CCC
%This just continues to loop through updating the parameters like above
%based on the resulting Jacobian and whether the new result is better than
%the previous step in parameter space
while SSE>1e3 && iter<maxiter && cccTC<1 
    disp(['Iteration ' num2str(iter)])
    if ~lambdaup
        %Recalculate the Jacobian every 25 iterations for efficiency
        if (~rem(iter,25)||recalc==1) 
            clear TCtest J* Del;
            TCtest = zeros((n^2)*pslices,numPs);
            J = zeros((n^2)*pslices,numPs);
            for dd = 1:numPs/comps
            parfor i = 1 + comps*(dd-1):comps + comps*(dd-1) 
                testprop = intermedparams;
                testprop(whichones(i)) = paramsnew(whichones(i));
                diff = testprop(whichones(i)) ...
                                            - intermedparams(whichones(i));
                if diff == 0
                    testprop(whichones(i)) = testprop(whichones(i))*0.99;
                    diff = testprop(whichones(i)) ...
                                            - intermedparams(whichones(i));
                end
                DMatrix = testprop(1)*ones(n,n,pslices);
                kMatrix = reshape(testprop(2:end),n,n,pslices);

                BigTCtest = NtcFDmech3DUTLogisticOnly(Z,imagedims,...
                            DMatrix,dt,kMatrix,CarryCaps,BCF,timesteps,...
                                                              G,DiffyG,nU);  
                TCtest(:,i) = reshape(BigTCtest,[],1);
                J(:,i) = (TCtest(:,i)-TCintermed)/diff;
            end
            end
            %Again just catching the remainder of parameters if any
            if mod(numPs,comps)~=0
            parfor i = 1+numPs-mod(numPs,comps):numPs                
                testprop = intermedparams;
                testprop(whichones(i)) = paramsnew(whichones(i));
                diff = testprop(whichones(i)) ...
                                            - intermedparams(whichones(i));
                if diff == 0
                    testprop(whichones(i)) = testprop(whichones(i))*0.99;
                    diff = testprop(whichones(i)) ...
                                            - intermedparams(whichones(i));
                end
                DMatrix = testprop(1)*ones(n,n,pslices);
                kMatrix = reshape(testprop(2:end),n,n,pslices);

                BigTCtest = NtcFDmech3DUTLogisticOnly(Z,imagedims,...
                            DMatrix,dt,kMatrix,CarryCaps,BCF,timesteps,...
                                                              G,DiffyG,nU);  
                TCtest(:,i) = reshape(BigTCtest,[],1);
                J(:,i) = (TCtest(:,i)-TCintermed)/diff;
            end
            end
        end
    end
    if(any(isnan(J)))
        disp('Jacobian is not an answer, ending optimization.'); 
        break 
    end
    
    Jt = J';
    Js = Jt*J;
    
    alphaval = (lam*trace(Js)*SSE^2);
    ncelldiffs = Zf-TCintermed;
    H = Js + sqrt(alphaval)*eye(size(Js,1),size(Js,2));

    g = (Jt)*ncelldiffs;
    clear Hsc gsc deltasc tempsc
    gsc = zeros(size(H,1),1); 
    Hsc = zeros(size(H,1),size(H,2));
    for i1 = 1:size(H,1)
        gsc(i1,1) = g(i1,1)/sqrt(H(i1,i1));
        for i2 = 1:size(H,2)
            Hsc(i1,i2) = H(i1,i2)/(sqrt(H(i1,i1))*sqrt(H(i2,i2))); 
        end
    end
    
    [L, U] = lu(Hsc);
    tempsc = L\gsc;
    delsc = U\tempsc;
    Del = zeros(size(delsc));
    for i1 = 1:size(delsc,1)
        Del(i1,1) = delsc(i1,1)/sqrt(H(i1,i1));
    end
    temp = zeros(size(intermedparams));
    temp(whichones) = Del;
    Del = temp;
    if(intermedparams(1)+Del(1)>uD)
        Del(1) = (uD-intermedparams(1))/2;
    end
    if(intermedparams(1)+Del(1)<0)
        Del(1) = -(intermedparams(1))/2;
    end
    
    newtest = intermedparams;
    newtest(whichones) = intermedparams(whichones) + Del(whichones);
    err = norm(Del)/norm(intermedparams(whichones));
    
    DMatrix = newtest(1)*ones(n,n,pslices);
    kMatrix = reshape(newtest(2:end),n,n,pslices);   

    BigTCnewtest = NtcFDmech3DUTLogisticOnly(Z,imagedims,DMatrix,dt,...
                              kMatrix,CarryCaps,BCF,timesteps,G,DiffyG,nU);    

    TCnewtest = reshape(BigTCnewtest,[],1); 
    SSEnew = sum((Zf-TCnewtest).^2);

    %If the new result is better than the previous we accept this as the
    %next step in parameter space, this is based on the SSE
    if SSEnew<SSE
        paramsnew = intermedparams;
        intermedparams = newtest;
        clear TCintermed;
        TCintermed = TCnewtest;
        cccTC = ccc_barnes2(Zf,TCnewtest);
        SSE = SSEnew; 
        if ~rem(iter,3) 
            lam = lam/9;
        end
        iter = iter + 1;
        lambdaup = 0;
        clear newtest SSEnew Jt Js H g gsc Hsc deltasc TCnewtest;
        upcnt = 0;
        recalc = 0;
    else %if we do not get a better result, update the lambda change
        if upcnt>10 && recalc==0
            lambdaup = 0;
            upcnt = 0;
            testcnt = testcnt + 1;
            recalc = 1;
            disp(['Recalculating Jacobian, recalculation count = ' ...
                                                        num2str(testcnt)]);
        else
            lam = lam*11;
            lambdaup = 1;
            upcnt = upcnt + 1;
            recalc = 0;
            disp(['Increasing lambda, lambda = ' num2str(lam)])
        end
    end
    if(testcnt>10||lam>10^50)  %fail safe
        disp('Lambda too big, ending optimization... '); 
        break 
    end
end
disp(['While loop final iteration ' num2str(iter)])

disp('Calculating final Jacobian')
TCtest = zeros((n^2)*pslices,numPs);
J = zeros((n^2)*pslices,numPs);
for dd = 1:numPs/comps
parfor i = 1 + comps*(dd-1):comps + comps*(dd-1) 
    testprop = intermedparams;
    testprop(whichones(i)) = paramsnew(whichones(i));
    diff = testprop(whichones(i)) - intermedparams(whichones(i));
    if diff == 0
        testprop(whichones(i)) = testprop(whichones(i))*0.99;
        diff = testprop(whichones(i)) - intermedparams(whichones(i));
    end
    DMatrix = testprop(1)*ones(n,n,pslices);
    kMatrix = reshape(testprop(2:end),n,n,pslices);

    BigTCtest = NtcFDmech3DUTLogisticOnly(Z,imagedims,DMatrix,dt,...
                              kMatrix,CarryCaps,BCF,timesteps,G,DiffyG,nU);  
    TCtest(:,i) = reshape(BigTCtest,[],1);
    J(:,i) = (TCtest(:,i)-TCintermed)/diff;
end
end
%Again just catching the remainder of parameters if any
if mod(numPs,comps)~=0
parfor i = 1+numPs-mod(numPs,comps):numPs
    testprop = intermedparams;
    testprop(whichones(i)) = paramsnew(whichones(i));
    diff = testprop(whichones(i)) - intermedparams(whichones(i));
    if diff == 0
        testprop(whichones(i)) = testprop(whichones(i))*0.99;
        diff = testprop(whichones(i)) - intermedparams(whichones(i));
    end
    DMatrix = testprop(1)*ones(n,n,pslices);
    kMatrix = reshape(testprop(2:end),n,n,pslices);

    BigTCtest = NtcFDmech3DUTLogisticOnly(Z,imagedims,DMatrix,dt,...
                              kMatrix,CarryCaps,BCF,timesteps,G,DiffyG,nU);  
    TCtest(:,i) = reshape(BigTCtest,[],1);
    J(:,i) = (TCtest(:,i)-TCintermed)/diff;
end
end

save(['params_' patientID '.txt'],'intermedparams','-ascii');
save(['TCfinal_' patientID '.txt'],'TCintermed','-ascii');
save(['err_' patientID '.txt'],'err','-ascii');
save(['SSE_' patientID '.txt'],'SSE','-ascii');
save(['ccc_' patientID '.txt'],'cccTC','-ascii');
save(['finaliteration_' patientID '.txt'],'iter','-ascii'); 
save(['Jacobian_' patientID '.txt'],'J','-ascii'); 


%end of file
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~