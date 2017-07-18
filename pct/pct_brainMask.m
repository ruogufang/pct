function mask = pct_brainMask(im,lb,ub,dsize)
% FUNCTION MASK = pct_brainMASK(IM, LB, UB, DSIZE) finds the brain mask on the given image
% by eliminating negative values and values exceeding the upper limit.
%
% INPUT:
%       IM      - Input image [Y x X]
%       LB      - Lower bound (default 0)
%       UB      - Upper bound (default 2000)
%       DSIZE   - Disk radius for morphological closing [Scalar] default:3
%
% OUTPUT:
%       MASK - Brain mask [Logical] (1 for brain region, 0 for non-brain
%       region)
%
% -Ruogu Fang 12/22/2010

if nargin < 2
    lb = 0;
    ub = 2000;
end

if nargin < 4
    dsize = 3;
end

cc = bwconncomp(lb < im & im <= ub);
numPixels = cellfun(@numel,cc.PixelIdxList);
[~,idx]=max(numPixels);
mask = false(size(im));
mask(cc.PixelIdxList{idx})=true;
% mask = false(size(im));
% mask(lb<im & im <=ub)=true;
str = strel('disk',dsize);
mask = imclose(mask,str);

mask = logical(mask);

end
    
