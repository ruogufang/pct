function pct_vol2gif(V,name)
%PCT_VOL2GIF converts 3D volumetric data to movies in gif format
%
%   This function generate GIF movie of the brain perfusion in cine mode from 3D volumetric data 
%
%   USAGE: PCT_VOL2GIF(V,NAME,SLICE);
%
%   PRE:
%       V       - A 3D volumetric brain perfusion data [T x Y x X]
%       NAME    - A name of output gif file [String]
%
%   Ruogu Fang 09/2011
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%

if nargin < 2 || isempty(name)
    name = 'CTP';
end

nframes = size(V,1);

I = squeeze(V(1,:,:));
c = [0 200];
ctshow(I,[],c,0);
axis tight
set(gca,'nextplot','replacechildren','visible','off')
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');

for k = 1 : nframes
  imshow(squeeze(V(k,:,:)),c);
  f = getframe;
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
end

im(:,:,1,nframes+1:nframes+2) = 0;

gif_name = [name '.gif'];
imwrite(im,map,gif_name,'DelayTime',1,'LoopCount',inf);
end
