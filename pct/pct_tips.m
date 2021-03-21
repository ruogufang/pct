% PCT_TIPS Two dimensional time-intensity profile similarity bilateral filtering.
%    This function implements 2-D TIPS bilateral filtering using the method
%    outlined in:
%
%    TIPS bilateral noise reduction in 4D CT perfusion scans produces
%    high-quality cerebral blood flow maps. Phys Med Biol. 2011 Jul
%    7;56(13):3857-72. doi: 10.1088/0031-9155/56/13/008. Epub 2011 Jun 8.
%
%    B = tips2(A,W,SIGMA) performs 2-D TIPS bilateral filtering. A should be a
%    double precision matrix of size NxMxT. The half-size of the Gaussian
%    bilateral filter window is defined by W. The standard deviations of the
%    bilateral filter are given by SIGMA, where the spatial-domain standard
%    deviation is given by SIGMA(1) and the intensity-domain standard deviation
%    is given by SIGMA(2).
%
%   Ruogu Fang 12/03/2014 Smadt Medical Informatics Learning and Evaluation
%   (SMILE) School of Computing and Information Sciences Florida International
%   University

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-process input and select appropriate filter.
function B = pct_tips(A,w,sigma)

% Verify that the input image exists and is valid.
if ~exist('A','var') || isempty(A)
   error('Input image A is undefined or invalid.');
end


% Verify bilateral filter window size.
if ~exist('w','var') || isempty(w) || ...
      numel(w) ~= 1 || w < 1
   w = 5;
end
w = ceil(w);

% Verify bilateral filter standard deviations.
if ~exist('sigma','var') || isempty(sigma) || ...
      numel(sigma) ~= 2 || sigma(1) <= 0 || sigma(2) <= 0
   sigma = [3 0.1];
end

% Apply TIPS bilateral filtering.
B = tipsflt(A,w,sigma(1),sigma(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements TIPS bilateral filter.
function B = tipsflt(A,w,sigma_d,sigma_r)

% Pre-compute Gaussian domain weights.
[X,Y] = meshgrid(-w:w,-w:w);
G = exp(-(X.^2+Y.^2)/(2*sigma_d^2));

% Rescale range variance (using maximum luminance).
% sigma_r = 100*sigma_r;

% % Create waitbar.
% h = waitbar(0,'Applying bilateral filter...');
% set(h,'Name','Bilateral Filter Progress');

% Apply bilateral filter.
dim = size(A);
B = zeros(dim);
for i = 1:dim(1)
   for j = 1:dim(2)
      
         % Extract local region.
         iMin = max(i-w,1);
         iMax = min(i+w,dim(1));
         jMin = max(j-w,1);
         jMax = min(j+w,dim(2));
         I = A(iMin:iMax,jMin:jMax,:);
      
         % Compute Gaussian range weights.
         d = I - repmat(A(i,j,:),[size(I,1) size(I,2) 1]);
         H = exp(-mean(d.^2,3)/(2*sigma_r^2));
      
         % Calculate bilateral filter response.
         F = H.*G((iMin:iMax)-i+w+1,(jMin:jMax)-j+w+1);
         norm_F = sum(F(:));
         B(i,j,:) = sum(sum(repmat(F,[1 1 size(I,3)]).*I))/norm_F;
                
   end
%    waitbar(i/dim(1));
end

% % Close waitbar.
% close(h);