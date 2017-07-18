function [ bbbmap, x, ymap, R ] = pct_bbbp(ctmap, cbvmap, AIF, dt, rho, first, last, mask)
%PCT_BBBP Computes a permeability map (BBBP).
%
%   This algorithm calculates blood-brain barrier (BBB) permeability
%   using the Patlak model.
%
%   USAGE: [BBBP X YMAP R] = PCT_BBBP(CTMAP,CBVMAP,AIF,DT,RHO,FIRST,LAST,MASK);
%
%   PRE:
%       CTMAP   - Untruncated, preprocessed CT sequence in HU [T x Y x X]
%       CBVMAP  - A CBV (Cerebral Blood Volume) map in mL/100g [Y x X]
%       AIF     - A preprocessed Arterial Input Function [T x 1]
%       DT      - Time interval between samples in seconds [Scalar]
%       RHO     - Average brain tissue density in g/mL [Scalar]
%       FIRST   - The first frame to include in the calculations [Scalar]
%       LAST    - The last frame to include in the calculations [Scalar]
%       MASK    - Optional parameter. A logical mask [Y x X] that indicates
%                 which pixels are to be processed (processes those that are
%                 TRUE).
%
%   POST:
%       BBBP    - A map of brain permeability in mL/100g/min [Y x X]
%       X       - The independent variable of the patlak plot [T x 1]
%       YMAP    - A map of dependent variables of the patlak plot [T x Y x X]
%       R       - A map of coefficients of determination (R^2)
%
%   Kolbeinn Karlsson 06/05/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

%Get mapsize
[len,height,width] = size(ctmap);

if nargin < 8
    mask = true(height, width);
end

if last == -1
    last = len;
end

    function [r2] = rsq(y,yfit)
        %Calculates the coefficient of determination (r-squared)
        yresid = y - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1) * var(y);
        r2 = 1 - SSresid/SStotal;
    end


%Pre-allocate the output variable
bbbmap = zeros([height width]);
R = zeros([height width]);
ymap = zeros(size(ctmap));


%Calculate the AIF integral for each time step
aif_sum = zeros(size(AIF));
for n = 1:length(AIF)
    aif_sum(n) = dt*trapz(AIF(1:n));
    %aif_sum(n) = sum(AIF(1:n));
end

%Preconditioning the AIF (We don't want to divide by zero)
AIF(AIF < 1) = 1;

%Calculate x, the independent variable of the Patlak plot
x = aif_sum ./ AIF;

%Calculate the y map, the dependent variables of the Patlak plot
for i = 1:height
    for j = 1:width
        %Pass through mask
        if mask(i,j)
            %Calculate the y-map
            RIF = squeeze(ctmap(:,i,j));
            ymap(:,i,j) = RIF ./ AIF * 1/rho - cbvmap(i,j)/100;
            %Calculate the permeability
            xx = x(first:last);
            yy = squeeze(ymap(first:last,i,j));
            bbbmap(i,j) = xx \ yy;
            %Calculate the coefficient of determination
            R(i,j) = rsq(yy,bbbmap(i,j)*xx);
        else
            %Set entries to zero
            ymap(:,i,j) = 0;
            R(i,j) = 0;
            bbbmap(i,j) = 0;
        end
    end
end

%Convert the permeability map from mL/g/s to mL/100g/min
bbbmap = bbbmap * 60 * 100;

%Filter negative values
bbbmap(bbbmap<0)=0;

% Remove nonbrain tissue
bbbmap(~mask) = 0;

end



