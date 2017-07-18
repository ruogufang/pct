function [X_den,P]=denoise_TTV(Xobs,lambda,l,u,P_init,pars)
%This function implements the FISTA method for 4D Tensor TV (TTV4D) denoising problems. 
%
% 2014-06-16 Add non-iso regularization parameters
%
% INPUT
% Xobs ..............................an observed noisy image (1-4D)
% lambda ............................ parameter [4x1 vector]
% pars.................................parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.epsilon ..................... tolerance for relative error used in
%                                                       the stopping criteria (Default=1e-4)
% pars.print ..........................  1 if a report on the iterations is
%                                                       given, 0 if the  report is silenced
% pars.tv .................................. type of total variation
%                                                      penatly.  'iso' for isotropic (default)
%                                                      and 'l1' for nonisotropic
%  
% OUTPUT
% X_den ........................... The solution of the problem 
%                                            min{||X-Xobs||^2+2*lambda*TV(X)}


%Define the Projection onto the box
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

% Assigning parameres according to pars and/or default values
flag=exist('pars', 'var');
if (flag&&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
% if (flag&&isfield(pars,'epsilon'))
%     epsilon=pars.epsilon;
% else
%     epsilon=1e-4;
% end
% if(flag&&isfield(pars,'print'))
%     prnt=pars.print;
% else
%     prnt=1;
% end
if(flag&&isfield(pars,'tv'))
    tv=pars.tv;
else
    tv='iso';
end

[m,n,p,q]=size(Xobs);

% clear P; clear R;
if(isempty(P_init))
    P{1}=zeros(m-1,n,p,q);    P{2}=zeros(m,n-1,p,q);    P{3}=zeros(m,n,p-1,q);  P{4}=zeros(m,n,p,q-1);
    R{1}=zeros(m-1,n,p,q);    R{2}=zeros(m,n-1,p,q);    R{3}=zeros(m,n,p-1,q);  R{4}=zeros(m,n,p,q-1);
else
    P{1}=P_init{1};    P{2}=P_init{2};      P{3}=P_init{3};  P{4}=P_init{4};
    R{1}=P_init{1};    R{2}=P_init{2};      R{3}=P_init{3};  R{4}=P_init{4};
end

tkp1=1;count=0;i=0;

D=zeros(m,n,p,q);%fval=inf;fun_all=[];
while((i<MAXITER)&&(count<5))
%    fold=fval;  
    i=i+1;    
%     Dold=D;    
    Pold=P;    
    tk=tkp1;
    D=project(Xobs-Lforward(R, m, n, p, q, lambda));
    Q=Ltrans(D, m, n, p, q);
    %%%%%%%%%%
    % Taking a step towards minus of the gradient
    P{1}=R{1}+1/(8*lambda(1)+eps)*Q{1};
    P{2}=R{2}+1/(8*lambda(2)+eps)*Q{2};
    P{3}=R{3}+1/(8*lambda(3)+eps)*Q{3};
    P{4}=R{4}+1/(8*lambda(4)+eps)*Q{4};
    
    %%%%%%%%%%
    % Peforming the projection step
    switch tv
        case 'iso'
            A=cat(1,P{1}.^2,zeros(1,n,p,q))+cat(2,P{2}.^2,zeros(m,1,p,q))+cat(3,P{3}.^2,zeros(m,n,1,q))+cat(4,P{4}.^2,zeros(m,n,p,1));
            A=sqrt(max(A,1));
            P{1}=P{1}./A(1:m-1,:,:,:); P{2}=P{2}./A(:,1:n-1,:,:); P{3}=P{3}./A(:,:,1:p-1,:); P{4}=P{4}./A(:,:,:,1:q-1);
        case 'l1'
            P{1}=P{1}./(max(abs(P{1}),1));
            P{2}=P{2}./(max(abs(P{2}),1));
            P{3}=P{3}./(max(abs(P{3}),1));
            P{4}=P{4}./(max(abs(P{4}),1));
        otherwise
            error('unknown type of total variation. should be iso or l1');
    end

    %%%%%%%%%%
    %Updating R and t
    tkp1=(1+sqrt(1+4*tk^2))/2;
    
    R{1}=P{1}+(tk-1)/(tkp1)*(P{1}-Pold{1});
    R{2}=P{2}+(tk-1)/tkp1*(P{2}-Pold{2});
    R{3}=P{3}+(tk-1)/tkp1*(P{3}-Pold{3});
    R{4}=P{4}+(tk-1)/tkp1*(P{4}-Pold{4});

    
%     re=norm(D-Dold,'fro')/norm(D,'fro');
%     if (re<epsilon)
%         count=count+1;
%     else
%         count=0;
%     end
    C=Xobs-Lforward(P, m, n, p, q, lambda);
    PC=project(C);
%     fval=-norm(C-PC,'fro')^2+norm(C,'fro')^2;
%     fun_all=[fun_all;fval];
end
X_den=D; iter=i;
return

function X=Lforward(P, m, n, p, q, lambda)

%       [m2,n2]=size(P{1});
%       [m1,n1]=size(P{2});
% 
%       if (n2~=n1+1)
%           error('dimensions are not consistent')
%       end
%       if(m1~=m2+1)
%           error('dimensions are not consistent')
%       end
% 
%       m=m2+1;
%       n=n2;

      X=zeros(m,n,p,q);
      X(1:m-1,:,:,:)=lambda(1)*P{1};
      X(:,1:n-1,:,:)=X(:,1:n-1,:,:)+lambda(2)*P{2};
      X(:,:,1:p-1,:)=X(:,:,1:p-1,:)+lambda(3)*P{3};
      X(:,:,:,1:q-1)=X(:,:,:,1:q-1)+lambda(4)*P{4};
      X(2:m,:,:,:)=X(2:m,:,:,:)-lambda(1)*P{1};
      X(:,2:n,:,:)=X(:,2:n,:,:)-lambda(2)*P{2};
      X(:,:,2:p,:)=X(:,:,2:p,:)-lambda(3)*P{3};
      X(:,:,:,2:q)=X(:,:,:,2:q)-lambda(4)*P{4};
   end
 
   function P=Ltrans(X, m, n, p, q)
%       [m,n]=size(X);
      P{1}=X(1:m-1,:,:,:)-X(2:m,:,:,:);
      P{2}=X(:,1:n-1,:,:)-X(:,2:n,:,:);
      P{3}=X(:,:,1:p-1,:)-X(:,:,2:p,:);
      P{4}=X(:,:,:,1:q-1)-X(:,:,:,2:q);
   end
end
