function [pm, rm] = pct_aath(pct,aif,dt,mask)
%PCT_ACTH Estimates IRF parameters of a PCT series using the ACTH model
%
%   USAGE: [PM, RM] = PCT_ACTH(PCT,AIF,DT,MASK);
%
%   INPUT:
%       PCT     - A PCT series [T x Y x X]
%       AIF     - Arterial Input Function [T x 1]
%       DT      - Sampling interval [1x1]
%       MASK    - A mask indicating which voxels to process [Y x X]
%
%   OUTPUT:
%       PM      - Parameter matrix. Contains the output. [4 x Y x X]
%       RM      - Residue matrix. Contains the fitting residue of each voxel.
%                 Not to be confused with the impulse residue function [Y x X].
%
%   This function computes the impulse residue function parameters for each
%   voxel in a PCT series. The PM will contain F,E,Ve, and Tc of each pixel.

%Get input size
[T H W] = size(pct);

%First find the Tc values we must generate.
Tcc = 0:dt:12;

%Optimization parameters
x0 = [0.05 0.1 0.1];        %Initial value
lb = 0.00001*[1 1 1]';      %Lower bounds
ub = [0.1 1 1]';            %Upper bounds

options = optimset('Algorithm','active-set','FinDiffType','forward',...
    'Display','notify','MaxSQPIter',10,'MaxIter',200,'MaxFunEvals',1000,...
    'FunValCheck','on','TolFun',1e-15,'TolX',1e-6);

%Length of Tcc
L = length(Tcc);

%Create result matrix
parameters = zeros(L,3); %Temporary variable for parameters
residue = zeros(L,1);    %Temporary variable for residue
pm = zeros(4,H,W); %Will hold the final result
rm = zeros(H,W);   %Will hold the residue

%Loop through every pixel
for i = 1:H
    for j = 1:W
        if mask(i,j)
            %Compute the parameters for the pixel
            for k = 1:L
                %objFunc = @(x) sqrt(sum((pct(:,i,j)-thfunc(dt,aif,x(1),...
                %   x(2),x(3), Tcc(k))).^2 )/length(aif));
                objFunc = @(x) pct_objfunc(pct(:,i,j),dt,aif,x(1),x(2), ...
                    x(3), Tcc(k));
                [parameters(k,:) residue(k)] = ...
                    fmincon(objFunc,x0,[],[],[],[],lb,ub,[],options);
            end
            %Find the best fit
            [C I] = min(residue);
            pm(:,i,j) = [parameters(I,:) Tcc(I)];
            rm(i,j) = residue(I);
        end
    end
end

end

