function output = TTV_FCSA(input,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% minimize reg*TV(y) + 0.5*||Ay-b||_2^2
%
%%% Sep. 16, 2009, Written by Junzhou Huang at Rutgers University
%%% Sep. 2013, Revised by Ruogu Fang to add an adaptive step size s
%%% 2014-6-16 Add non-iso regularization parameter lambda [1x3]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parsin.MAXITER=1; parsin.tv='l1'; l=input.l; u=input.u; 

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
n1=input.n1; n2=input.n2; nt=input.nt;
A=input.A; b=input.b;
TV_iso = false; TV_eps = 0;
y=zeros(input.nt,input.n1*input.n2);
f=inf;
yr=zeros(size(y));
tnew=1;
reg = input.reg;
if isfield(input,'dt')
    dt = input.dt;
else
    dt = 1;
end

if isfield(input,'tt')
tt=input.tt;
else tt = nt;
end

% regularization parameter should be a 1x3 vector
switch length(reg)
    case 1
        reg = ones(1,3)*reg;
    case 2
        reg = [reg(1) ones(1,2)*reg(2)];
end

% iterno=1;
% output.funv(iterno)=0.5*sum(b(:).*b(:));
% output.xrate(iterno)=1;
% output.snr(iterno)=0;
% output.xtime(iterno)=0;

iterno=0;
s = zeros(size(y));

for itr = 1:input.maxitr  % total iter counter
    iterno=iterno+1;
    %     iitr = iitr+1;
    told=tnew;
    yp=y;
    %     yg=yr-(A'*(A*yr)-Atb)/Lx;
    parsin.s = s;
    
    % Steepest gradient descent
    e = A'*(b-A*yr);
    ev = e(:);
    Ae = A*e;
    Ae = Ae(:);
    step = (ev'*ev)/(Ae'*Ae);
    yg = yr + e*step;
    
    yg = reshape(yg,nt,n1,n2);
    yt = yg(1:tt,:,:);
    
    if (itr==1)
        [yt, P]=denoise_TTV(yt, 2*reg,-inf,inf,[],parsin);
    else
        [yt, P]=denoise_TTV(yt, 2*reg,-inf,inf,P,parsin);
    end
    
    y = yg; y(1:tt,:,:) = yt;
    y = reshape(y,nt,[]);
    y = project(y);
    grad_y = [Grad1(y) Grad2(y) Grad3(y)];
    sqsum_y = sum(grad_y.*grad_y,2);
    %     s = (reg^2*sqsum_y)/(reg^2*sqsum_y+reg2^2);
    
    [f, f_prev]=get_f(A*y, y, f);
    
    if f>f_prev
        y=yp;
        f=f_prev;
    end
    
    tnew=(1+sqrt(1+4*told^2))/2;
    yr = y+((told-1)/tnew)*(y-yp);
    if itr>30, tnew=1; end
    
    %     output.rel(iterno)=norm(y-yp, 'fro')/norm(yp, 'fro');
    output.funv(iterno)=f;
    %     output.snr(iterno)=snr(y, input.f);
    output.xtime(iterno)=cputime-t00;
    %     if iterno>2 if output.rel(iterno)>output.rel(iterno-1); tnew=1; end
    %     end
    
    if itr>=input.no
        break
    end
end

output.y=y/dt;

return;

    function [f, f_prev]=get_f(Ay, y, f)
        v1 = TV(y, reg);
        r3 = Ay - b; v3 = 0.5*(norm(r3).^2);
        f_prev = f;
        f = v1 + v3;
    end

    function p=Grad1(u)
        % backward finite difference along dim 1
        u = reshape(u,nt,n1,n2);
        p = cat(1,u(1,:,:),diff(u,1,1));
        p = p(:);
    end

    function q=Grad2(u)
        % backward finite difference along dim 2
        u = reshape(u,nt,n1,n2);
        q = cat(2,u(:,1,:),diff(u,1,2));
        q = q(:);
    end

    function q=Grad3(u)
        % backward finite difference along dim 2
        u = reshape(u,nt,n1,n2);
        q = cat(3,u(:,:,1),diff(u,1,3));
        q = q(:);
    end

    function TVx=TV(x, reg)
        % Total Variation norm of x, x is a n by n matrix
        grad_x = [reg(1)*Grad1(x) reg(2)*Grad2(x) reg(3)*Grad3(x)];
        if TV_iso
            pt_sqsum = sum(grad_x.*grad_x,2);
            if TV_eps == 0; TVx = sum(sqrt(pt_sqsum)); else TVx = sum(sqrt(pt_sqsum+TV_eps)); end
        else
            TVx = norm(grad_x(:),1);
        end
    end
end