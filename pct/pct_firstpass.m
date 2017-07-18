function [ bbbmap ] = pct_firstpass(ctmap, cbvmap, ttp, AIF, rho, first, last)
%PCT_FIRSTPASS Returns a BBB map with delay correction
%
%   USAGE: BBBMAP = PCT_FIRSTPASS(CTMAP,CBV,TTP,AIF,RHO,FIRST,LAST);
%
%   PRE:
%       CTMAP   - A pre-processed CT sequence [T x Y x X]
%       CBV     - A CBV map [Y x X]
%       TTP     - A TTP map [Y x X]
%       AIF     - An arterial input function [T x 1]
%       RHO     - Average brain tissue density in g/mL [Scalar]
%       FIRST   - The first frame to include in the calculations [Scalar]
%       LAST    - The last frame to include in the calculations [Scalar]
%
%   POST:
%       BBBMAP  - A map of blood-brain barrier permeability [Y x X]
%
%   Kolbeinn Karlsson 06/12/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%Start timer
tikky = tic;

%Truncate the data if necessary
if nargin == 6
   ctmap = ctmap(first:last,:,:);
   AIF = AIF(first:last);
end

%Pre-allocate the output variable
[height width] = size(cbvmap);
bbbmap = zeros([height width]);

%Calculate the AIF integral for each time step
aif_sum = zeros(size(AIF));
for n = 1:length(AIF)
    aif_sum(n) = sum(AIF(1:n));
end

%Calculate the TTP of the AIF
[V AIFTTP] = max(AIF);

%Calculate the BBB permeability
for j = 1:height
    for i = 1:width
        AIFDC = pct_timedelay(AIF, ttp(i,j)-AIFTTP);
        X = aif_sum ./ AIFDC;
        RIF = squeeze(ctmap(:,j,i));
        Y = RIF ./ AIFDC * 100/rho - cbvmap(j,i);
        bbbmap(j,i) = X \ Y;
    end
end

%Convert the permeability map from mL/100g/s to mL/100g/min
bbbmap = bbbmap .* 60;

%Show elapsed time
toc(tikky)

end


