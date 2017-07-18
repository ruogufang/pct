function [out] = pct_mask(in, mask)
%PCT_MASK Segments out tissue not in the mask
%   
%   Ruogu Fang 07/11/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  OUT = PCT_MASK(IN, MASK)
%   
%   PRE:
%       IN     - a CT series in Hounsfield Units [T x X x Y]
%       MASK   - Mask of brain tissue [Y x X]
%
%
%   POST:
%       OUT    - A CT map like IN but with voxels with false value in the
%       mask segmented out (changed to zero)
%
%   All voxels with value false in the mask are set to zero.

mask = permute(mask,[3 1 2]);
out = zeros(size(in));

for t = 1 : size(in,1)
    frame = squeeze(in(t,:,:));
    frame(~mask) = 0;
    out(t,:,:) = frame;
end

end

