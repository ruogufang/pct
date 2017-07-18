function [ mask ] = pct_segment2(img, loth, hith )
%PCT_SEGMENT2 Segments an image using thresholds and region growing
%
%   USAGE:  MASK = PCT_SEGMENT2(IMG, LOTH, HITH);
%
%   INPUT:
%       IMG     - A 2D input image [Y x X]
%       LOTH    - The lower threshold [Scalar]
%       HITH    - The higher threshold [Scalar]
%
%   OUTPUT:
%       MASK    - A logical mask defining a central, connected region within the
%                 defined thresholds [Y x X]
%
%   This function segments an image using thresholding and returns only the
%   center-most connected region within the thresholding bounds.
%
%   Kolbeinn Karlsson, 08/14/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

stack=java.util.Stack();
siz = size(img);
mask = false(siz);

%Threshold the image
tmask = img > loth & img < hith;

%Locate a central region
x = round(siz(1)/2);
y = round(siz(2)/2);
while tmask(x,y) == 0
    y = y+1;
end

%Now grow the region
stack.push([x y]);

while ~stack.isEmpty
    p = stack.pop;
    if mask(p(1),p(2))
        continue
    end
    mask(p(1),p(2)) = 1;
    if p(1) > 1 && tmask(p(1)-1,p(2))
        stack.push([p(1)-1,p(2)]);
    end   
    if p(1) < siz(1) && tmask(p(1)+1,p(2))
        stack.push([p(1)+1,p(2)]);
    end
    
    if p(2) > 1 && tmask(p(1),p(2)-1)
        stack.push([p(1),p(2)-1]);
    end   
    if p(2) < siz(2) && tmask(p(1),p(2)+1)
        stack.push([p(1),p(2)+1]);
    end

end

