function [AIF,aif_x,aif_y] = pct_aifsearch(data,r,range,text)
%PCT_AIFSEARCH finds the arterial input function from CTP data
%
%   Ruogu Fang Revised 08/22/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  AIF = PCT_AIFSEARCH(DATA,R);
%
%   PRE:
%       data    - CTP input data [T x Y x X]
%       r       - Radius for AIF or VOF searching around pivot point [Scalar]
%       range   - display range [lo hi] (optional)
%       text    - text displayed in the graphical window [String] (optional)
%
%   POST:
%       AIF     - Arterial input function [T x 1]
%       AIF_X   - Column coordinate of AIF [Scalar]
%       AIF_Y   - Row coordinate of AIF [Scalar]
%
% PCT_AIF finds the optimal AIF or VOF within the searching radius r of
% pivot point designated by the user. The optimal AIF or VOF is the one
% with the highest peak value on the attenuation curve.
%

if nargin < 3
    range = [0 50];
end

if nargin < 4
    text = 'Please select the artery';
end

%Get size of CTP data
[time,height,width] = size(data);

%Extract the middle frame in the time series
I = squeeze(data(round(time/2),:,:));

%Get pivot point for AIF searching
figure;imshow(I,range); title(text);
[px,py] = ginput(1);
px = round(px);
py = round(py);

% Draw white rectangle on the searching area
hold on;
rectangle('Position',[px-r,py-r,2*r,2*r],'LineWidth',3,'EdgeColor','r');

%All candidate attenuation curves within radius r of the pivot point
curves = data(:,py-r:py+r,px-r:px+r);

%Find the peaks of attenuation curves
peaks = squeeze(max(curves,[],1));

%Find the optimal AIF index with the highest peak
[maxPeak,aifIndex] = max(peaks(:));

%Convert index to subscripts
[aif_y, aif_x] = ind2sub(size(peaks),aifIndex);

%Find the optimal AIF
AIF = curves(:,aif_y,aif_x);

%Show the result curve
figure;plot(AIF);

%Ask user if the curve is good to continue the computation
choice = questdlg('Would you like to continue with this curve?',...
    'CTP Deconvolution');
% Handle response
switch choice
    case 'Yes'
        close all;
    case 'No'
        close;
        clf;
        AIF = pct_aifsearch(data,r,range);
    case 'Cancel'
        close all;
end

aif_x = aif_x + px - r - 1;
aif_y = aif_y + py - r - 1;

end


