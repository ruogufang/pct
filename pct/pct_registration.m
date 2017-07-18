function outmap = pct_registration(ctmap)
%PCT_REGISTRATION Registers the CT images to the first frame
%
%   This algorithm registers the CT images using rigid transformation.
%
%   USAGE: [OUTMAP] = PCT_REGISTRATION(CTMAP);
%
%   PRE:
%       CTMAP   - Untruncated, unpreprocessed CT sequence in HU [T x Y x X]
%
%   POST:
%       OUTMAP  - Reigstered CT maps
%
%   Ruogu Fang 01/30/13
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%

[len height width] = size(ctmap);
outmap = zeros(size(ctmap));
fixed = squeeze(ctmap(1,:,:));
transformType = 'rigid';
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 300;
optimizer.MinimumStepLength = 5e-4;
outmap(1,:,:) = ctmap(1,:,:);

for i = 2:len
    im = squeeze(ctmap(i,:,:));
    imreg = imregister(im,fixed,transformType,optimizer,metric);
    outmap(i,:,:) = imreg;
end

end