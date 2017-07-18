function bb = pct_minBoundingBox(I,square)
%PCT_MINBOUNDINGBOX finds the minimum bounding box with boundaries along
%the axis
%
%   USAGE: BB = PCT_MINBOUNDINGBOX(I);
%
%   PRE:
%       I       - A logical image [Y x X]
%       squre   - Flag indicating output a square [0 or 1]
%
%   POST:
%       BB      - A matrix containing the coordinates of the bounding box
%                   corners [2 x 4]
%
%   Ruogu Fang 10/05/13
%   Advanced Multimedia Processing (AMP) Lab, Cornell University

if nargin < 2
    square = 0;
end

[r,c]=find(I);

bb = zeros(2,4);

if square == 0
    bb(:,1) = [min(r);min(c)];
    bb(:,2) = [max(r);min(c)];
    bb(:,3) = [max(r);max(c)];
    bb(:,4) = [min(r);max(c)];
else
    w = max(c) - min(c);
    h = max(r) - min(r);
    l = round(max(w,h)/2);
    m = round([1/2*(max(r)+min(r)),1/2*(max(c)+min(c))]);
    bb(:,1) = [m(1)-l;m(2)-l];
    bb(:,2) = [m(1)+l;m(2)-l];
    bb(:,3) = [m(1)+l;m(2)+l];
    bb(:,4) = [m(1)-l;m(2)+l];
end

end