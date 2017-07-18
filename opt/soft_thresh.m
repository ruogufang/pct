% Function to solve soft thresholding problem
%
% arg min_{x} ||x - b||_{2}^{2} + lambda*||x||_{1}
%
% Usage:- x = soft_thresh_cvx
%
% where:- <in>
%         b = bias vector
%         lambda = weighting on the l1 penalty
%         <out>
%         x = solution          
%
% Written by Simon Lucey 2012
% or
% soft_thresh = @(b,lambda) sign(b).*max(abs(b) - lambda/2,0);

function x = soft_thresh(b,lambda)

% Set the threshold
th = lambda/2; 

% First find elements that are larger than the threshold
k = find(b > th);
x(k) = b(k) - th; 

% Next find elements that are less than abs
k = find(abs(b) <= th); 
x(k) = 0; 

% Finally find elements that are less than -th
k = find(b < -th); 
x(k) = b(k) + th; 
x = x(:); 