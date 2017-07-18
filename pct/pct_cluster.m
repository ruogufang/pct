function [tissues masks] = pct_cluster(data,cbf,thresh)
% Cluster PCT data into different tissue types based on temporal average
% enhancement
%
% Input: 
%   data -  CTP data [TxYxX]
%   cbf - CBF map [YxX]
%   thresh - segmentation thresholds (vector [(N+1)x1]) (N= # tissue types)
%
% Output:
%   tissues - cell array of cbf maps for each tissue type [Nx1] cell
%   masks - cell array of logical masks for each tissue type [Nx1] cell
%
% Ruogu Fang
% 9/22/2012

N = length(thresh)-1;
meanmap = squeeze(mean(data));

tissues = cell(N,1);
masks = cell(N,1);

types = length(thresh)-1;
cluster = zeros(size(cbf));
cluster(cbf>thresh(2)) = cbf(cbf>thresh(2));
tissues{1} = cluster;
mask = zeros(size(cbf));
mask(cbf>thresh(2)) = 1;
masks{1} = mask;
cluster = zeros(size(cbf));
cluster(cbf<=thresh(2) & meanmap>thresh(3)) = cbf(cbf<=thresh(2) & meanmap>thresh(3));
tissues{2} = cluster;
mask = zeros(size(cbf));
mask(cbf<=thresh(2) & meanmap>thresh(3)) = 1;
masks{2} = mask;

for i = 3 : types
    cluster = zeros(size(cbf));
    mask = zeros(size(cbf));
    cluster(meanmap <= thresh(i) & meanmap > thresh(i+1) & cbf<=thresh(2)) = cbf(meanmap<=thresh(i) & meanmap>thresh(i+1) & cbf<=thresh(2));
    mask(meanmap <= thresh(i) & meanmap > thresh(i+1) & cbf<=thresh(2)) = 1;
    tissues{i} = cluster;
    masks{i} = mask;
%     figure;imshow(cluster);
end

end

