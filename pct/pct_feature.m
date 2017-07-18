function feature = pct_feature(data,mask)
% Dynamic feature of PCT data
%
% Pre:
%   data    - 3D PCT volume [T x X x Y]
%   mask    - binary mask where 1 indicates data and 0 indicate nonbrain [X x Y]
%
% Post:
%   feature - feature map extracted from data [X x Y x 3]
%
% Ruogu Fang 2/22/2013

%Median enhancement (ME)
ME = pct_windowing(squeeze(median(data)),0,100);
ME(~mask) = 0;

%Peak enhancement (PE)
PE = pct_windowing(squeeze(max(data)),0,150);
PE(~mask) = 0;
 
%Area under curve (AUC)
AUC = pct_windowing(squeeze(sum(data)),0,10000);
AUC(~mask) = 0;

feature = cat(3,ME,PE,AUC);
