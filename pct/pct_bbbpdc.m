function [bbbmap,xmap,ymap,R] = pct_bbbpdc(ctmap,cbvmap,ttpmap,AIF,dt,rho,first,last,mask)
%PCT_BBBPDC Computes a delay-corrected permeability map (BBBP).
%correction
%
%   This algorithm calculates blood-brain barrier permeability (BBBP)
%   using the Patlak model and applying the Wintermark delay correction.
%
%   USAGE: [BBBP XMAP YMAP R] = PCT_BBBPDC(CTMAP,CBVMAP,TTPMAP,AIF,DTRHO,FIRST,LAST);
%
%   PRE:
%       CTMAP   - Untruncated, preprocessed CT sequence in HU [T x Y x X]
%       CBVMAP  - A CBV (Cerebral Blood Volume) map in mL/100g [Y x X]
%       TTPMAP  - A map of Time-To-Peak in seconds [Y x X]
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
%       XMAP    - A map of independent variables of the patlak plot [T x 1]
%       YMAP    - A map of dependent variables of the patlak plot [T x Y x X]
%       R       - A map of coefficients of determination (R^2) [Y x X]
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

    function [r2] = rsq(y,yfit)
        %Calculates the coefficient of determination (r-squared)
        yresid = y - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1) * var(y);
        r2 = 1 - SSresid/SStotal;
    end

    function [y] = integrate(x,dt)
        %Returns the integral (sum) of x where dt is the interval between
        %samples
        y = zeros(size(x));
        for p = 1:length(x)
            %Use right Riemann sums for fast computation
            %             y(p) = dt*sum(x(1:p));
            %Use Trapz numerical integration
            y(p) = dt*trapz(x(1:p));
        end
    end


%Pre-allocate the output variable
[height width] = size(cbvmap);
bbbmap = zeros([height width]);
R = zeros([height width]);
ymap = zeros(size(ctmap));
xmap = zeros(size(ctmap));


%Preconditioning the AIF (we don't want to divide by zero)
AIF(AIF <= 1) = 1;

%Get the TTP(AIF)
ttp_aif = find(AIF == max(AIF));

AIF_sum_orig = integrate(AIF,dt);

for i = 1:height
    for j = 1:width
        if mask(i,j)
            %Get the delayed AIF and the integral of AIF
            delay = ttpmap(i,j)*1/dt-ttp_aif;
            if delay < 0
                delay = 0;
            end
            AIF_temp = pct_timedelay(AIF,delay);
            AIF_sum = pct_timedelay(AIF_sum_orig,delay);
            %Calculate the x variable
            xmap(:,i,j) = AIF_sum ./ AIF_temp;
            %Calculate the y variable
            RIF = squeeze(ctmap(:,i,j));
            ymap(:,i,j) = RIF ./ AIF_temp * 1/rho - cbvmap(i,j)/100;
            %Calculate the permeability
            xx = squeeze(xmap(first+delay:last,i,j));
            yy = squeeze(ymap(first+delay:last,i,j));
            bbbmap(i,j) = xx \ yy;
%             if bbbmap(i,j) > 0.001
%                 plot(xx,yy);                
%             end
            %Calculate the coefficient of determination
            R(i,j) = rsq(yy,bbbmap(i,j)*xx);
        end
    end
end

%Convert the permeability map from mL/g/s to mL/100g/min
bbbmap = bbbmap * 60 * 100;

end



