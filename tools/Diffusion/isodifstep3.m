function y = isodifstep(x, d)
% ISODIFSTEP   Isotropic diffusion step
%
%    y = ISODIFSTEP(x, d) calculates de isotropic (scalar) diffusion
%    step "y" based on the image "x" and on the diffusivity "d". If
%    "d" is constant the diffusion will be linear, if "d" is
%    a matrix the same size as "x" the diffusion will be nonlinear.
%
%    The diffused image is calculated as:
%      xd = x + T*isodifstep(x,d)  , T = step size
%

% Translations of d
%dop = roll(d,[0 1]);
%dom = roll(d,[0 -1]);
%dpo = roll(d,[1 0]);
%dmo = roll(d,[-1 0]);
d_dpoo = d + roll(d,[ 1   0   0]);
d_dmoo = roll(d_dpoo,[-1   0   0]);
d_dopo = d + roll(d,[ 0   1   0]);
d_domo = roll(d_dopo,[ 0  -1   0]);
d_doop = d + roll(d,[ 0   0   1]);
d_doom = roll(d_doop,[ 0   0  -1]);


% Translations of x
%xop = roll(x,[0 1]);
%xom = roll(x,[0 -1]);
%xpo = roll(x,[1 0]);
%xmo = roll(x,[-1 0]);
x_xpoo = x - roll(x,[ 1   0   0]);
x_xmoo = roll(x_xpoo,[-1   0   0]); % Must multiply -1
x_xopo = x - roll(x,[ 0   1   0]);
x_xomo = roll(x_xopo,[ 0  -1   0]); % Must multiply -1
x_xoop = x - roll(x,[ 0   0   1]);
x_xoom = roll(x_xoop,[ 0   0  -1]); % Must multiply -1


% Calculate y = dx/dt
%y = -.5 * ( (d+dmo).*(x-xmo) + (dpo+d).*(x-xpo) + (d+dom).*(x-xom) + (dop+d).*(x-xop)  );
%y =   .5 * ( (d_dmo).*(x_xmo) - (d_dpo).*(x_xpo) + (d_dom).*(x_xom) - (d_dop).*(x_xop)  );
y =   .5 * ( (d_dmoo).*(x_xmoo) - (d_dpoo).*(x_xpoo) + (d_domo).*(x_xomo) - (d_dopo).*(x_xopo)  + (d_doom).* (x_xoom)- (d_doop).*(x_xoop));