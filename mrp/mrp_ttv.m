function [y,funv] = mrp_ttv(input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MRP_TTV minimize reg*TV(x) + 0.5*||Ax-b||_2^2
%
%   USAGE:  [Y,FUNV] = MRP_TTV(INPUT,VARARGIN)
%
%   INPUT.
%       x       - Spatio-temporal residue functions [T x N]
%       A       - Teopliz matrix of arterial input function [T x T]
%       b       - Contrast concentration curves [T x N]
%       reg     - anisotropic regularization parameter [4 x 1]
%       maxiter - Maximum number of iterations
%       no      - Number of iterations
%       tt      - Time length for TV regularization [Scalar]
%       dim     - Dimension of the actual x [ndim x 1]
%       tv      - TV type: iso or l1
%       l       - Lower bound
%       u       - Upper bound
%       L       - Parameter for TV denoise
%       MASK    - A logical mask. The computation will only be performed
%                 for the voxels that are logical 1. [T x X x Y x Z]
%
%   OUTPUT:
%       Y       - BF * R. The impulse residue function scaled by the cerebral
%                 blood flow. [T x X X Y X Z].
%       FUNV    - Cost function over iterations.
%
%   This function solves the TV regularized Indicator-Dilution equation
%
%       C = F * conv( C_art,R ) + reg * TV(R)
%   
%
%%% Sep. 16, 2009, Written by Junzhou Huang at Rutgers University
%%% Sep. 2013, Revised by Ruogu Fang to add an adaptive step size s
%%% 2014-6-16 Add non-iso regularization parameter lambda [1x3]
% 2014-11-13 Extend to 4D, Ruogu Fang
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parsin.MAXITER=1; l=input.l; u=input.u; 

if((l==-Inf)&&(u==Inf))
    project=@(x)x;
elseif (isfinite(l)&&(u==Inf))
    project=@(x)(((l<x).*x)+(l*(x<=l)));
elseif (isfinite(u)&&(l==-Inf))
    project=@(x)(((x<u).*x)+((x>=u)*u));
elseif ((isfinite(u)&&isfinite(l))&&(l<u))
    project=@(x)(((l<x)&(x<u)).*x)+((x>=u)*u)+(l*(x<=l));
else
    error('lower and upper bound l,u should satisfy l<u');
end
t00 = cputime;

A=input.A; b=input.b; dim = input.dim;

if ~isfield(input,'reg')
    reg = ones(length(dim),1);
else
    reg = input.reg; 
end
if ~isfield(input,'mask')
    mask = ones(dim);
else
    mask = input.mask;
end
if ~isfield(input,'tv')
    tv = 'l1';
else
    tv = input.tv;
end
parsin.tv=tv; 
    
y=zeros([dim(1) prod(dim(2:end))]);
f=inf;
yr=zeros([dim(1) prod(dim(2:end))]);
tnew=1;

if isfield(input,'dt')
    b = b/input.dt; % correct for sampling rate
end

if isfield(input,'tt')
    tt=input.tt;
else
    tt = dim(1);
end

% regularization parameter should be a 1x3 vector
switch length(reg)
    case 1
        reg = ones(1,4)*reg;
    case 2
        reg = [reg(1) ones(1,4)*reg(2)];
end

iterno=0;
s = zeros(size(y));

for itr = 1:input.maxitr  % total iter counter
    iterno=iterno+1;
    told=tnew;
    yp=y;
    parsin.s = s;
    
    % Steepest gradient descent
    e = A'*(b-A*yr);
    ev = e(:);
    Ae = A*e;
    Ae = Ae(:);
    step = (ev'*ev)/(Ae'*Ae);
    yg = yr + e*step;
    
    yg = reshape(yg,dim);
    yt = yg(1:tt,:,:);
    
    if (itr==1)
        [yt, P]=denoise_TTV(yt, 2*reg,-inf,inf,[],parsin);
    else
        [yt, P]=denoise_TTV(yt, 2*reg,-inf,inf,P,parsin);
    end
    
    y = yg; y(1:tt,:,:) = yt;
    y = reshape(y,dim(1),[]);
    y = project(y);
    %     grad_y = [Grad1(y) Grad2(y) Grad3(y)];
    %     sqsum_y = sum(grad_y.*grad_y,2);
    %     s = (reg^2*sqsum_y)/(reg^2*sqsum_y+reg2^2);
    
    [f, f_prev]=get_f(A*y, y, f, b, reg, tv);
    
    if f>f_prev
        y=yp;
        f=f_prev;
    end
    
    tnew=(1+sqrt(1+4*told^2))/2;
    yr = y+((told-1)/tnew)*(y-yp);
    if itr>30, tnew=1; end
    
    %     output.rel(iterno)=norm(y-yp, 'fro')/norm(yp, 'fro');
    funv(iterno)=f;
    %     output.snr(iterno)=snr(y, input.f);
    xtime(iterno)=cputime-t00;
    %     if iterno>2 if output.rel(iterno)>output.rel(iterno-1); tnew=1; end
    %     end
    
    if itr>=input.no
        break
    end
    
end

y = reshape(y,dim);
mask = shiftdim(repmat(mask,[ones(length(dim)-1,1);dim(1)]),length(dim)-1);
y(~mask) = 0;

return;

    function [f, f_prev]=get_f(Ay, y, f, b, reg, tv)
        v1 = TV(y,reg,tv);
        r3 = Ay - b; v3 = 0.5*(norm(r3).^2);
        f_prev = f;
        f = v1 + v3;
    end

    function grad_x = Grad(u,reg)
        % backward finite difference in all dimensions
        grad_x = [];
        for d = 1 : ndims(u)
            u1 = shiftdim(u,d-1);
            u1 = u1(1,:,:,:);
            u1 = shiftdim(u1,ndims(u)-d+1);
            p = cat(d,u1,diff(u,1,d));
            grad_x = [grad_x reg(d)*p(:)];
        end
    end

    function TVx = TV(x,reg,tv)
        % Total Variation norm of x, x is a high-dimensional tensor
        grad_x = Grad(x,reg);
        switch tv
            case 'iso'
                pt_sqsum = sum(grad_x.*grad_x,2);
                TVx = sum(sqrt(pt_sqsum));
            case 'l1'
                TVx = norm(grad_x(:),1);
            otherwise
                error('unknown type of total variation. should be iso or l1');
        end
    end

end