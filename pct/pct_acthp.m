function [out ] = pct_acthp(Q,AIF,dt)
%function [F,E,Ve,Tc ] = pct_acthp(Q,AIF,dt)
%PCT_ACTHP Deconvolves the I-D eqn using the ACTH model for one pixel
%
%   USAGE: [F,E,VE,TC] = PCT_ACTHP(Q,AIF,DT);
%
%   INPUT: 
%       Q       - A Tissue-attenuation curve [T x 1]
%       AIF     - Arterial Input Function [T x 1]
%       DT      - The sampling interval [1 x 1]
%
%   OUTPUT:
%       F       - Blood flow [1x1]
%       E       - Extraction fraction [1x1]
%       VE      - Extravascular blood volume [1x1]
%       TC      - Transit time through capillary [1x1]
%
%   This function devonvolves the indicator-dilution equation for a single pixel
%   using the Adibatic Correction to the Tissue Homogeneity model. It uses the
%   SQP algorithm to find the best fit parameters for different discrete values
%   of TC, then returns the set of values with the best fit. The TC values tried
%   are the ones between 0-16.
%
%   Kolbeinn Karlsson, 09/14/12
%   Advanced Multimedia Processing (AMP) Lab, Cornell University



%First find the Tc values we must generate. 
Tcc = 0:dt:12;

%Optimization parameters
x0 = [0.05 0.1 0.1];        %Initial value
lb = [0.00001 0.00001 0.00001]';      %Lower bounds
ub = [0.1 1 1]';              %Upper bounds

% options = optimset('Algorithm','active-set','TolX',0.00000001,...
%      'FinDiffType','central', 'Display','notify','TolFun',0.00000001,...
%      'MaxFunEvals',1000); 
 
options = optimset('Algorithm','active-set','FinDiffType','forward',...
    'Display','notify','MaxSQPIter',10,'MaxIter',100,'MaxFunEvals',1000,...
    'FunValCheck','on','TolFun',1e-15,'TolX',1e-6); 

%Length of Tcc
L = length(Tcc);
 
%Create result matrix
parameters = zeros(L,3);
residue = zeros(L,1);
 
%Now run SQP for different values of Tcc
for i = 1:L
    objFunc = @(x) sqrt(sum((Q-thfunc(dt,AIF,x(1),x(2),x(3), ...
        Tcc(i))).^2 )/length(AIF));
    [parameters(i,:) residue(i)] = ...
        fmincon(objFunc,x0,[],[],[],[],lb,ub,[],options);
end

%Find the best fit
[C I] = min(residue);

% %Return the parameters
% F = parameters(I,1);
% E = parameters(I,2);
% Ve = parameters(I,3);
% Tc = Tcc(I);

out.para = parameters;
out.res = residue;
out.tc = Tcc';
out.ind = I;








end

