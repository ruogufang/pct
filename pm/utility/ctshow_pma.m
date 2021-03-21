%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [c,cm,cb] = ctshow_pma(im, mask, c, cflag, CLT)
% displays CT image by using CLT (PMA colormap) and rescaling the region in the mask
%
% OUTPUTS:
%   c       - Output gray value range to be displayed
%   cm      - Output colormap
%   cb      - Output colorbar
%
% INPUTS:
%   im      - Images to be scaled and displayed
%   mask    - region of interest in the image to be scaled and displayed
%   c       - User specified gray value range to be displayed
%   cflag   - Indicator for graylevel/color display
%               default:default colormap, pma:PMA colormap
%   CLT     - Color Lookup Table
%
% 2/14/2019 - Yao Xiao @ SMILE | UF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c,cm,cb] = ctshow_pma(im,mask,c,cflag,CLT)

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
    cflag = 'default';
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

switch cflag
    case 'default'
        colormap('default');
        cmap = colormap;
        cmap(1,:) = [0 0 0];
        cm = colormap(cmap);
        cb = colorbar;
    
    case 'pma'
        cm = colormap(CLT);
        cb = colorbar;
end

end

