function out = mrp_gamfit(in,mask)
%MRP_GAMFIT fits the concentration-time curve (CTC) with Gamma-variant function
%
%   Ruogu Fang 11/05/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University
%
%   USAGE:  OUT = MRP_GAMFIT(IN, MASK);
%
%   PRE:
%       IN     - Concentration-time curves [T x X x Y x Z]
%       MASK   - Brain Mask [X x Y x Z]
%
%   POST:
%       OUT    - Gamma-variant fitted Concentration-time curves [T x X x Y x Z]
%
%   This function fits the concentration-time curves with Gamma-variant
%   functions in the form of
%   v(t) = c0 * (t-t0)^alpha*e^(-(t-t0)/beta)
%
% TO BE COMPELETED...

[T,X,Y,Z] = size(in);

if nargin < 2 || isempty(mask)
    loth = 0;
    hith = 1;
    PRE = 2:10;
    B = squeeze(mean(in(PRE,:,:,:)));
    mask = pct_brainMask(B,loth,hith);
end

out = zeros(size(in));
p0 = [2,0,3,1.8];

gamvarFun =  @(p,x) (p(1) .* (x-p(2)).^ p(3).*exp(-(x-p(2))./(p(4)+eps)));

for i = 1 : X
    for j = 1 : Y
        for k = 1 : Z
            if mask(i,j,k);
                C = squeeze(in(:,i,j,k));
                C(C<=0) = 1;
                p = nlinfit(1:T,C',gamvarFun,p0);
                Cg = gamvarFun(p,1:T);
                out(:,i,j,k) = Cg;
            end
        end
    end
end

end


