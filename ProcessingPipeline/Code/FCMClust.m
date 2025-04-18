function [center, U, obj_fcn] = FCMClust(data, cluster_n, options)  
% FCMClust.m   ????C??????data??cluster_n?  
%  
% ???  
%   1.  [center,U,obj_fcn] = FCMClust(Data,N_cluster,options);  
%   2.  [center,U,obj_fcn] = FCMClust(Data,N_cluster);  
%    
% ???  
%   data        ---- nxm??,??n???,??????m?????  
%   N_cluster   ---- ??,????????,????  
%   options     ---- 4x1?????  
%       options(1):  ?????U????>1                  (???: 2.0)  
%       options(2):  ??????                           (???: 100)  
%       options(3):  ????????,??????           (???: 1e-5)  
%       options(4):  ????????????                (???: 1)  
% ???  
%   center      ---- ????  
%   U           ---- ?????  
%   obj_fcn     ---- ?????  
%   Example:  
%       data = rand(100,2);  
%       [center,U,obj_fcn] = FCMClust(data,2);  
%       plot(data(:,1), data(:,2),'o');  
%       hold on;  
%       maxU = max(U);  
%       index1 = find(U(1,:) == maxU);  
%       index2 = find(U(2,:) == maxU);  
%       line(data(index1,1),data(index1,2),'marker','*','color','g');  
%       line(data(index2,1),data(index2,2),'marker','*','color','r');  
%       plot([center([1 2],1)],[center([1 2],2)],'*','color','k')  
%       hold off;  
  
   
  
if nargin ~= 2 & nargin ~= 3,    %???????????2??3?  
 error('Too many or too few input arguments!');  
end  
  
data_n = size(data, 1); % ??data????(rows)?,?????  
in_n = size(data, 2);   % ??data????(columns)????????  
% ??????  
default_options = [2; % ?????U???  
    100;                % ??????  
    1e-5;               % ????????,??????  
    1];                 % ????????????  
  
if nargin == 2,  
 options = default_options;  
 else       %???options????????  
 % ??????????????????option;  
 if length(options) < 4, %??????opition???4?????????;  
  tmp = default_options;  
  tmp(1:length(options)) = options;  
  options = tmp;  
    end  
    % ??options??????0(?NaN),?????1  
 nan_index = find(isnan(options)==1);  
    %?denfault_options???????????options???????.  
 options(nan_index) = default_options(nan_index);  
 if options(1) <= 1, %?????????????1  
  error('The exponent should be greater than 1!');  
 end  
end  
%?options ?????????????;  
expo = options(1);          % ?????U???  
max_iter = options(2);  % ??????  
min_impro = options(3);  % ????????,??????  
display = options(4);  % ????????????  
  
obj_fcn = zeros(max_iter, 1); % ???????obj_fcn  
  
U = initfcm(cluster_n, data_n);     % ?????????,?U???????1,  
% Main loop  ????  
for i = 1:max_iter,  
    %??k??????????ceneter,?????U?????;  
 [U, center, obj_fcn(i)] = stepfcm(data, U, cluster_n, expo);  
 if display,  
  fprintf('FCM:Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));  
 end  
 % ??????  
 if i > 1,  
  if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro,  
            break;  
        end,  
 end  
end  
  
iter_n = i; % ??????  
obj_fcn(iter_n+1:max_iter) = [];  
  
  
% ???  
function U = initfcm(cluster_n, data_n)  
% ???fcm????????  
% ??:  
%   cluster_n   ---- ??????  
%   data_n      ---- ????  
% ???  
%   U           ---- ?????????  
U = rand(cluster_n, data_n);  
col_sum = sum(U);  
U = U./col_sum(ones(cluster_n, 1), :);  
  
   
  
% ???  
function [U_new, center, obj_fcn] = stepfcm(data, U, cluster_n, expo)  
% ??C??????????  
% ???  
%   data        ---- nxm??,??n???,??????m?????  
%   U           ---- ?????  
%   cluster_n   ---- ??,????????,????  
%   expo        ---- ?????U???                       
% ???  
%   U_new       ---- ?????????????  
%   center      ---- ????????????  
%   obj_fcn     ---- ?????  
mf = U.^expo;       % ?????????????  
center = mf*data./((ones(size(data, 2), 1)*sum(mf'))'); % ?????(5.4)?  
dist = distfcm(center, data);       % ??????  
obj_fcn = sum(sum((dist.^2).*mf));  % ??????? (5.1)?  
tmp = dist.^(-2/(expo-1));      
U_new = tmp./(ones(cluster_n, 1)*sum(tmp));  % ????????? (5.3)?  
  
   
  
% ???  
function out = distfcm(center, data)  
% ??????????????  
% ???  
%   center     ---- ????  
%   data       ---- ???  
% ???  
%   out        ---- ??  
out = zeros(size(center, 1), size(data, 1));  
for k = 1:size(center, 1), % ????????  
    % ??????????????????????  
    out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));  
end  