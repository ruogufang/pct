function [out,S0] = mrp_convert(in,mask,TE)
%MRP_CONVERT Converts MR signal values to contrast concentration
%
%   Ruogu Fang 11/05/2014
%   Smart Medical Informatics Learning and Evaluation (SMILE)
%   School of Computing and Information Sciences
%   Florida International University
%
%   USAGE:  OUT = MRP_CONCONV(IN, FRAMES);
%
%   PRE:
%       IN     - MRP data [T x X x Y x Z]
%       MASK   - Brain Mask [X x Y x Z]
%       TE     - Echo time [Sclar] (sec) default: 0.04 sec
%
%   POST:
%       OUT    - A series of contrast concentration maps [T x X x Y x Z]
%
%   This function calcualtes the average of each pixel in the first
%   FRAMES frames as the precontrast signal intensity and use logorithm relation
%   to convert the MR signal intensity to the contrast concentration
%
% The conversion relationship is:
% C (concentration of contrast) = -1/(k*TE*rho) * ln(S/S0)

k = 1; % experiment specific constant assumed common for all tissues
PRE = 2:10;
    
dim = ndims(in);
[T,X,Y,Z] = size(in);

if nargin < 3
    TE = 40; % echo time of the acquisition sequence
end

if nargin < 4
    rho = 1.05; % density of brain tissue
end

if nargin < 2 || isempty(mask)
    loth = 100;
    hith = 4500;
    B = squeeze(mean(in(PRE,:,:,:)));
    mask = pct_brainMask(B,loth,hith);
end

% Calculate the precontrast signal intensity
S0 = median(in(PRE,:,:,:));

% in = in(:,:,:,delay+1:end);
out = zeros(size(in));
S0 = repmat(S0,[T 1]);
mask_v = shiftdim(repmat(mask,[ones(1,dim-1) T]),dim-1);

out(mask_v) = -k/TE .* log(in(mask_v)./(S0(mask_v)+eps));

out(isinf(out))=eps;
out(out<0) = eps;

end


