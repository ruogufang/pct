function [AIF,aif_x,aif_y] = pct_aifautoselect(data,mask,roi)
%PCT_AIFAUTOSELECT finds the arterial input function from CTP data automatically
%
%   Ruogu Fang 10/16/2014
%   School of Computing and Information Sciences
%   Florida International University
%
%   USAGE:  AIF = PCT_AIFAUTOSELECT(DATA,ROI);
%
%   PRE:
%       data    - CTP input data [T x Y x X]
%       roi     - Region of interest for AIF selection [4x1 vector (x_1,y_1,x_2,y_2)]
%
%   POST:
%       AIF     - Arterial input function [T x 1]
%       AIF_X   - Column coordinate of AIF [Scalar]
%       AIF_Y   - Row coordinate of AIF [Scalar]
%
% PCT_AIFAUTOSELECT finds the optimal AIF within the searching radius r of
% pivot point designated by the user. The optimal AIF or VOF is the one
% with the highest peak value on the attenuation curve.
%

% If no mask input, mask is whole image
if nargin < 2
    mask = ones(size(data,2),size(data,3));
end

% If no ROI input, ROI is the entire 2D image size
if nargin < 3
    roi = [1,1,size(data,2),size(data,3)];
end

T_cutoff = 60; % consider only T_cutoff time frames due to recirculation
thresh = 0.2; % percentage of maximum peak value to be considered for AIF

%Get size of CTP data
[time,height,width] = size(data);
T = min(T_cutoff,time);

%Extract the middle frame in the time series
I = squeeze(data(round(time/2),:,:));

%Get pivot point for AIF searching
figure;imshow(I,[0 50]);

% Draw white rectangle on the searching area
hold on;
rectangle('Position',[roi(1),roi(2),roi(3)-roi(1),roi(4)-roi(2)],'LineWidth',3,'EdgeColor','r');

%All candidate attenuation curves within radius r of the pivot point
curves = data(1:T,roi(1):roi(3),roi(2):roi(4));

%Find the peaks of all attenuation curves
[peaks,tmax] = max(curves,[],1);
peaks = squeeze(peaks);
tmax = squeeze(tmax);

%Sort peak values and index
[sortPeaks,idxPeaks] = sort(peaks(:),'descend');

%Remove peak values less than a threshold
idx = find(sortPeaks<sortPeaks(1)*thresh & ~mask(:));
sortPeaks(idx) = [];
idxPeaks(idx) = [];
tmax_thresh = tmax(idxPeaks);

%Sort time to peak values for the highest peaks
[sortTmax, idxTmax] = sort(tmax_thresh(:),'ascend');

%Convert index to subscripts
[aif_y, aif_x] = ind2sub(size(peaks),idxTmax(1));

%Find the optimal AIF
AIF = curves(:,aif_y,aif_x);

aif_x = aif_x + roi(1) - 1;
aif_y = aif_y + roi(2) - 1;

% Draw red point at AIF location found in the brain slice
plot(aif_x,aif_y,'o','LineWidth',2);

%Show the result curve
figure;plot(AIF);



end


