function [ out ] = bilateralct(img, slice_no, n, sigma_s, sigma_t, sw, mask)
%BILATERALCT Filters a CT series using a 4D bilateral filter
%
%   USAGE:  OUT = BILATERALCT(IMG, Z, SLICE_NO, SIGMA_S, SIGMA_T, SW, MASK);
%
%   PRE:
%       IMG     - A 4D array [Y x X x Z x T]
%       SLICE_NO- The desired slice [Scalar]
%       N       - Filter size. Must be all odd numbers [M, P, Q, R]
%       SIGMA_S - 2-element array containing the domain std (1) and the range
%                 std (2) for the spatial filtering [1 x 2]
%       SIGMA_T - 2-element array containing the domain std (1) and the range
%                 std (2) for the temporal filtering [1 x 2]
%       SW      - (optional parameter) The dimensional weights [A x B x C x D].
%                 Recommended value for CT: [1 1 2 1].
%       MASK    - (optional parameter) A logical mask the same size as the first
%                 three dimensions of IMG [Y x X x Z]. The filter will only
%                 process the voxels that are logical true.
%
%   POST:
%       OUT     - The output slice [Y x X x Z]
%
%   This function is intended for use with 4D (3D + time) CT data. It takes a 4D
%   data set as input, but outputs only a single slice (2D + time). A 4D
%   bilateral filter with Gaussian weighting functions is applied to the slice.
%   The 4D bilateral filter is approximated by applying first a 3D bilateral
%   filter on the spatial dimensions followed by a 1D bilateral filter on the
%   temporal dimension.
%
%   Kolbeinn Karlsson, 07/25/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University


%Get the image size
[H W L T] = size(img);

%Precalculate constants with respect to the loop
S_sr = 2*sigma_s(2)^2;
S_tr = 2*sigma_t(2)^2;

%%Create Gaussian weighting function for spatial filtering
%Get the half-sizes
hs_x = (n(1)-1)/2;
hs_y = (n(2)-1)/2;
hs_z = (n(3)-1)/2;
hs_t = (n(4)-1)/2;
%Create the meshgrid
[x y z] = meshgrid(-hs_x:hs_x,-hs_y:hs_y,-hs_z:hs_z);
%Impart spatial weights on the grid
x = x*sw(1);
y = y*sw(2);
z = z*sw(3);
%Create the Gaussian
gaussian_sd = exp(-(x.^2 + y.^2 + z.^2)/(2*sigma_s(1)^2));
%Apply the bounds for gaussian_sd
kmin = max(1,slice_no-hs_z);
kmax = min(L,slice_no+hs_z);
Gd = gaussian_sd(:,:,(kmin:kmax)-slice_no+hs_z+1);

%Pre-allocate output image
out = zeros(H,W,T);

%%Apply the spatial filter
for t = 1:T
    for i = 1+hs_x:H-hs_x
        for j = 1+hs_y:W-hs_y
            if mask(i,j,slice_no)
                %Get the 3D window function
                I = img(i-hs_x:i+hs_x,j-hs_y:j+hs_y,kmin:kmax,t);
                %Compute range weighting function
                Gr = exp(- (I-img(i,j,slice_no,t)).^2/S_sr );
                %Apply the spatial filter
                Ws = Gr .* Gd;
                out(i,j,t) = sum(Ws(:).*I(:))/sum(Ws(:));
            end
        end
    end
end

%%Create Gaussian weighting function for the temporal filtering
t_index = -hs_t:hs_t;
%Apply the weight
t_index = sw(4)*t_index;
%Create the gaussian
gaussian_td = exp(- (t_index.^2)/(2*sigma_t(1)^2));

for i = 1+hs_x:H-hs_x
    for j = 1+hs_y:W-hs_y
        if mask(i,j,slice_no)
            for t = 1:T
                %Get the local window
                tmin = max(1,t-hs_t);
                tmax = min(T,t+hs_t);
                I = out(i,j,tmin:tmax);
                %Create the range Gaussian
                Gr = exp(- (I-out(i,j,t)).^2/S_tr);
                %Truncate the domain Gaussian if necessary
                Gd = gaussian_td( (tmin:tmax)-t+hs_t+1 );
                Wi = Gr(:) .* Gd(:);
                out(i,j,t) = sum(Wi .* I(:)) / sum(Wi);
            end
        end
    end
end
                


end

