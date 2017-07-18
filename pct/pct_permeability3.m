function [ bbbmap, xmap, ymap, rmse ] = pct_permeability3(ctmap, cbvmap, ttpmap, AIF, rho, first, last, mask)
%PCT_PERMEABILITY2 Calculates a permeability map from a CBV map with delay
%correction
%
%   This algorithm calculates blood-brain barrier (BBB) permeability
%   using the Patlak model and applying the delay correction.
%
%   USAGE: [BBBP XMAP YMAP R] = PCT_PERMEABILITY2(CTMAP,CBVMAP,TTPMAP,AIF,RHO,FIRST,LAST);
%
%   PRE:
%       CTMAP   - Untruncated, preprocessed CT sequence in HU [T x Y x X]
%       CBVMAP  - A CBV (Cerebral Blood Volume) map in mL/100g [Y x X]
%       TTPMAP  - A map of Time-To-Peak in seconds [Y x X]
%       AIF     - A preprocessed Arterial Input Function [T x 1]
%       RHO     - Average brain tissue density in g/mL [Scalar]
%       FIRST   - The first frame to include in the calculations [Scalar]
%       LAST    - The last frame to include in the calculations [Scalar]
%       MASK    - Optional parameter. A logical mask [Y x X] that indicates
%                 which pixels are to be processed (processes those that are
%                 TRUE).
%
%   POST:
%       BBBP    - A map of brain permeability in mL/100g/min [Y x X]
%       XMAP    - A map of independent variables of the patlak plot [T x 1]
%       YMAP    - A map of dependent variables of the patlak plot [T x Y x X]
%       RMSE    - A map of root mean square errors (RMSE) [Y x X]
%
%   Kolbeinn Karlsson 06/05/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%

%Get map dimensions
[height width] = size(cbvmap);

%If mask is not supplied, have no mask
if nargin < 8
    mask = true(height, width);
end

    function [m] = rootmse(y,yfit)
        %Calculates the root mean square error
        yresid = y - yfit;
        SSresid = sum(yresid.^2);
        m = sqrt(SSresid/length(yresid))/(max(y)-min(y));
    end

    function [y] = integrate(x)
        %Returns the integral (sum) of x
        y = zeros(size(x));
        for kuk = 1:length(x) %I use the strange name to avoid naming conflicts
            y(kuk) = sum(x(1:kuk));
        end
    end
        

%Pre-allocate the output variable
[height width] = size(cbvmap);
bbbmap = zeros([height width]);
rmse = zeros([height width]);
ymap = zeros(size(ctmap));
xmap = zeros(size(ctmap));

%Calculate the AIF integral for each time step
%aif_sum = zeros(size(AIF));
%for n = 1:length(AIF)
%    aif_sum(n) = sum(AIF(1:n));
%end

%Preconditioning the AIF (we don't want to divide by zero)
AIF(AIF <= 1) = 1;

%Get the TTP(AIF)
ttp_aif = find(AIF == max(AIF));

%Precondition the TTP
% lo = round(0.75*ttp_aif);
% hi = round(1.5*ttp_aif);
% ttpmap(ttpmap < lo) = ttp_aif;
% ttpmap(ttpmap > hi) = ttp_aif;
%ttpmap(ttpmap < lo) = lo;
%ttpmap(ttpmap > hi) = hi;




for i = 1:height
    for j = 1:width
        if mask(i,j)
            %Get the delayed AIF and the integral of AIF
            AIF_temp = pct_timedelay(AIF,ttpmap(i,j)-ttp_aif);
            AIF_sum = integrate(AIF_temp);
            %Calculate the x variable
            xmap(:,i,j) = AIF_sum ./ AIF_temp;
            %Calculate the y variable
            RIF = squeeze(ctmap(:,i,j));
            ymap(:,i,j) = RIF ./ AIF_temp * 100/rho - cbvmap(i,j);
            %Calculate the permeability
            xx = squeeze(xmap(first:last,i,j));
            yy = squeeze(ymap(first:last,i,j));
            bbbmap(i,j) = xx \ yy;
            %Calculate the coefficient of determination
            rmse(i,j) = rootmse(yy,bbbmap(i,j)*xx);
        end
    end
end

%Convert the permeability map from mL/100g/s to mL/100g/min
bbbmap = bbbmap .* 60;

end



