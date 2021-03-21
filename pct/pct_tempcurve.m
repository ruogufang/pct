function [x y c] = pct_tempcurve(data)
%PCT_TEMPCURVE plots the temporal tissue enhancement curve at selected
%pixel
%
%   Ruogu Fang 07/11/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  [X Y] = PCT_TEMPCURVE(DATA);
%
%   PRE:
%       data    - CTP input data [Y x X x 1 x T]
%
%   POST:
%       x       - x coordinate of the select pixel
%       y       - y coordinate of the select pixel
%       c       - temporal TEC curve of the selected pixel
%
%   This function requires user to select a point in the displayed image
%   and will show the temporal curve of the specificed pixel. If user is
%   not satisfied with the selection, left click on the plot; otherwise
%   right click on the plot and the function ends.
%

data = squeeze(permute(data,[4 1 2 3]));
I = double(squeeze(data(1,:,:)));
button = 1;
while button == 1
    imshow(I,[0 100]);
    [x y] = ginput(1);
    
    y = round(y);
    x = round(x);
    
    c = data(:,y,x);
    plot(c);
    
    [x1 y1 button] = ginput(1);
end
fprintf('x=%d,y=%d\n',x,y);
end

