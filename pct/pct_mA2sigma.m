function sigma = pct_mA2sigma(mA, mA0, K)
%PCT_MA2SIGMA Computes the corresponding sigma value for a specific mA
%
%   Ruogu Fang 7/7/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  SIGMA = PCT_MA2SIGMA(MA);
%
%   PRE:
%       MA        - Simulated mA value of tube current in human CT
%       [Scalar]
%       MA0       - True mA value of tube current in human CT
%       [Scalar]
%       K         - Parameter controlling the relationship between sigma
%                   and inverse of square root of mA (default K=103.09 for
%                   human brain)
%
%   POST:
%       SIGMA     - Standard deviation of Gaussian noise [Scalar]
%
%   This function calculates the standard variation of Gaussian noise at a
%   simulated tube current value in mA
%

if nargin < 3
    K = 103.09;
end

sigma = K * sqrt(1/mA - 1/mA0);

end