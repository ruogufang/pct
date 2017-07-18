function tissues = pct_cluster_cbf(CBF,thresh)
% Cluster PCT data into different tissue types
% Ruogu Fang
% 9/22/2012

tissues = cell(0);

types = length(thresh)-1;
for i = 1 : types
    cluster = zeros(size(CBF));
    cluster(CBF < thresh(i) & CBF > thresh(i+1)) = CBF(CBF<thresh(i) & CBF>thresh(i+1));
    tissues{i} = cluster;
end

end