function [ bbbpmap ] = pct_bbbp(rifmap, cbvmap, aif, rho, pre, post)
%PCT_BBBP Calculates a permeability map (BBBP).
%
%   Ruogu Fang  06/19/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  BBBPMAP = PCT_BBBP(RMAP, RHO)
%
%   PRE:
%       rifmap  - A map of region input functions [T x Y x X]
%       cbvmap  - A map of cerebral blood volume [Y x X]
%       aif     - Artery input function [T x 1]
%       rho     - Brain tissue density in g/mL [Scalar]
%       pre     - First frame to calculate BBBP [Scalar]
%       post    - Last frame to calculate BBBP [Scalar]
%
%   POST:
%       bbbpmap  - A map of blood-brain barrier permeability in
%       mL/100g/min [Y x X]
%
%   This function calculates the blood-brain barrier permeability from a
%   map of residue functions and arterial input functions using Patlak
%   model. The BBBP is calculated as the slope in Y = CBV + K * X, where Y
%   = RIF(t)/AIF(t) and X = int(AIF(tau)d_tau)/AIF(t) using linear
%   regression.

[T Y X] = size(rifmap);
aif = smooth(aif,'moving');
[max_aif,idx_aif] = max(aif);
bbbpmap = zeros(Y,X);

cbvmap = cbvmap * rho / 100;

for i = 1 : Y
    for j = 1 : X
        if cbvmap(i,j) ~= 0
            rif = smooth(rifmap(:,i,j));
            [max_rif,delay] = max(rif(idx_aif:60));            
            civ = zeros(T,1);
            if delay > 0
                civ(delay:end) = aif(1:end-delay+1);
            else
                civ = aif;
            end
            x = sum(civ) ./ (civ+eps);
            y = rifmap(:,i,j) ./ (civ+eps);
            pre1 = max(pre,pre+delay);
            x = x(pre1:post);
            y = y(pre1:post);
            p = x\(y-cbvmap(i,j));
            bbbpmap(i,j) = p(1);
        end
    end
end

scaling_factor_bbbp = 60 * 100 / rho;
bbbpmap = scaling_factor_bbbp * bbbpmap;

end

