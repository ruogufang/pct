function [px,py,regs] = pct_lcurve(A,b,sz,regs)

%PCT_LCURVE generates the L-curve for parameter selection of reg 
% for the problem
% x = argmin (1/2||Ax-b||_2^2+reg*||x||_TV)
%

if nargin < 4
    regs = [0.1:0.2:1 2:10]; % ranage of candidate regs
end
N = length(regs);
px = zeros(N,1);
py = zeros(N,1);

% FCSA parameters
input.n1=sz(2);input.n2=sz(3); input.nt = sz(1);
input.maxitr=500;
input.num=1; 
input.l=-inf; input.u=inf;
input.no = 50;
input.A = A;
input.b = b;

for i = 1 : N
    reg = regs(i);
    input.reg = reg;
    out = TVR_FCSA(input);
    x = out.y;
    px(i) = norm(A*x-b,'fro')/(sz(2)*sz(3));
    r = reshape(x,sz);
    py(i) = pct_tv(r)/prod(sz);
    
end

end
