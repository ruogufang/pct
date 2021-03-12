function mA = pct_sigma2mA(sigma)
%PCT_SIGMA2MA Computes the corresponding mA value for a specific sigma
%
%   Ruogu Fang 06/21/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  MA = PCT_SIGMA2MA(SIGMA);
%
%   PRE:
%       SIGMA     - Standard deviation of Gaussian noise [Scalar]
%
%   POST:
%       MA        - Corresponding mA value of tube current in human CT
%       [Scalar]
%
%   This function calculates the simulated tube current value in mA given
%   the standard deviation (sigma) of the added Gaussian noise
%

mA = 190*103.09^2/(190*sigma^2+103.09^2);

end