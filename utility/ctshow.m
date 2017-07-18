%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function c = ctshow(im, mask, c, cflag)
% displays CT image by rescaling the region in the mask
%
% INPUTS:
%   im      - Images to be scaled and displayed
%   mask    - region of interest in the image to be scaled and displayed
%   c       - User specified gray value range to be displayed
%   cflag   - Indicator for graylevel/color display
%               0(default):gray, 1:color
%
%9/21/2010 - Ruogu Fang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = ctshow(im,mask,c,cflag)

%Read inputs
if nargin < 2 || isempty(mask)
    mask = true(size(im,1),size(im,2));
end
if nargin < 3 || isempty(c)
    thresh = 0.05;
    [n,xout] = hist(double(im(mask&im>0)),1000);
    r = find(n >= max(n) * thresh);
    c = real([xout(r(1)) max(xout(r(end)),xout(r(1)+1))]);
end
if nargin < 4
    cflag = 1;
end

%Show image
im(~mask) = 0;
if ismatrix(im)
    imshow(im,c);
else
    im = permute(im,[1 2 4 3]);
    montage(im,'DisplayRange',c);
end
axis tight;
axis equal;

if cflag
    colormap('default');
    cmap = colormap;
    cmap(1,:) = [0 0 0];
    colormap(cmap);
%     cbar_axes = colorbar;
%         set(cbar_axes,'YColor',[1 1 1]);
end

% impixelinfo;
