function y = anidifstep(x, a, b, c)
% ANIDIFSTEP   Anisotropic diffusion step
%
%    y = ANIDIFSTEP(x, a, b, c) calculates de anisotropic (tensor) 
%    diffusion step "y" based on the image "x" and on the diffusion tensor "D".
%       D = / a  b \
%           \ b  c /
%
%    The diffused image is calculated as:
%      xd = x + T*anidifstep(x,d)  , T = step size
%

% Translations of a
a_apo = a + roll(a,[1 0]);
a_amo = roll(a_apo,[-1 0]);
%a_apo = a + roll(a,[1  0]);
%a_amo = roll(a_apo,[-1 0]);

% Translations of b
bop = roll(b,[0 1]);
bom = roll(b,[0 -1]);
bpo = roll(b,[1 0]);
bmo = roll(b,[-1 0]);


% Translations of c
%cop = roll(c,[0  1]);
%com = roll(c,[0 -1]);
c_cop = c + roll(c,[0 1]);
c_com = roll(c_cop,[0 -1]);

% Translations of x
xop = roll(x,[ 0  1]);
xom = roll(x,[ 0 -1]);
xpo = roll(x,[ 1  0]);
xmo = roll(x,[-1  0]);
xpp = roll(x,[ 1  1]);
xpm = roll(x,[ 1 -1]);
xmp = roll(x,[-1  1]);
xmm = roll(x,[-1 -1]);

% Calculate y = dx/dt
%y = -.25*(bmo+bop).*xmp + .5*(cop+c).*xop + .25*(bpo+bop).*xpp + .5*(amo+a).*xmo - .5*(amo+2*a+apo+com+2*c+cop).*x ...
%    + .5*(apo+a).*xpo + .25*(bmo+bom).*xmm + .5*(com+c).*xom - .25*(bpo+bom).*xpm;

y = .5* ( (c_cop).*xop + (a_amo).*xmo - (a_amo + a_apo + c_com + c_cop).*x + (a_apo).*xpo + (c_com).*xom) ...
   + .25* ( -1*( (bmo+bop).*xmp + (bpo+bom).*xpm ) + (bpo+bop).*xpp + (bmo+bom).*xmm );