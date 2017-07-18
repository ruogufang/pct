function [ p, ci, r ] = pct_lincon(x, y, mask)
%PCT_LINCON computes Lin's concorndance correlation coefficient with 95%
%confidence interval
%
%   Ruogu Fang 04/17/2013 Advanced Multimedia Processing (AMP) Lab, Cornell
%   University
%
%   USAGE:  [p, CI] = PCT_LINCON(X, Y);
%
%   PRE:
%       X, Y    - Input vectors for correlation computation [Nx1]
%       MASK    - Valid numbers [N x 1 logical]
%
%   POST:
%       P       - Lin's concordance correlation coefficient [Scalar] 
%       CI      - 95% confidence interval [2x1 vector]
%       R       - Pearson's Correlation Coefficient [Scalar]
%
% Note: When N<30, bootstrap 2-sided 95% confidence bounds are calculated.
% To be added.
%
% See: 
%
% 1. Lawrence I-Kuei Lin (March 1989). "A concordance correlation
% coefficient to evaluate reproducibility". Biometrics (International
% Biometric Society) 45 (1): 255-268. 
% 2. Lawrence I-Kuei Lin (March 2000). "A Note on the Concordance
% Correlation Coefficient". Biometrics 56: 324?325.
%

if length(x) ~= length(y)
    fprintf('Warning: Input equal length vectors\n');
    p = NaN; ci = NaN;
    return
end

if nargin < 3
    mask = true(size(x));
end

x = x(mask);
y = y(mask);

N = length(x);
% if N < 30
%     fprintf('Warning: For Lin Concordance correaltion coefficient: N>=30\n');
%     p = NaN; ci = NaN;
%     return
% end

p = (2*(x-mean(x))'*(y-mean(y))/N)/(var(x,1)+var(y,1)+(mean(x)-mean(y))^2); % Lin'e concorndance correlation coefficient
r = ((x-mean(x))'*(y-mean(y))/N)/(std(x,1)*std(y,1)); % Pearson product moment correlation
se = sqrt(((1-r^2)/r^2*p^2*(1-p^2)+2*p^3*(1-p)*((mean(x)-mean(y))^2)/...
    (std(x,1)*std(y,1)*r)-p^4*(mean(x)-mean(y))^4/(2*std(x,1)^2*std(y,1)^2*r^2))/(N-2)); % standard error
ci = zeros(2,1);
ci(1) = p-1.96*se/sqrt(N);
ci(2) = p+1.96*se/sqrt(N);

end

