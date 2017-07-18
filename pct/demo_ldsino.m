% demo_ldsino
% test code of pct_ldsino.m
% simulate low-dose sinogram
%
% Reference:
% Zeng, D., Huang, J., Bian, Z., Niu, S., Zhang, H., Feng, Q., Liang, Z. and Ma, J., 2015. A simple low-dose x-ray CT simulation from high-dose scan. IEEE transactions on nuclear science, 62(5), pp.2226-2233.

% Ruogu Fang
% 8/6/2016

close all; clear; clc;

% User-defined phantom
E = [1200 0.69 0.92 0 0 0;
    -480 0.6224 0.874 0 -0.0184 0;
    -120 0.41 0.16 -0.22 0 108;
    -120 0.31 0.11 0.22 0 72;
    60 0.25 0.21 0 0.35 90;
    60 0.046 0.046 0 0.1 0;
    60 0.023 0.046 0 -0.1 0;
    80 0.046 0.023 -0.08 -0.605 0;
    80 0.023 0.023 0 -0.605 0;
    80 0.046 0.023 0.06 -0.605 90;
    36 0.029 0.029 0 0.61 0];

% create phantom
xtrue = phantom(512);
% xtrue = xtrue/max(xtrue(:));
xtrue = rot90(xtrue);

f.down = 1;
f.do_sparse = 1; % set to one for sparse angle case
ig = image_geom('nx', 512, 'fov', 50, 'down', f.down);
ig.mask = ig.circ > 0;
sg = sino_geom('ge1', 'units', 'cm', 'strip_width', 'd', ...
    'down', f.down);
if f.do_sparse, sg.na = 50; end % # views in sparse angle case

% generate sinogram from phantom data
A = Gtomo2_strip(sg,ig);
% A = Gtomo2_dscmex(sg, ig);
hsino = A*xtrue;

% show high-dose image and sinogram
im plc 2 2
clim = [0 0.4];
im(1, xtrue, 'x', clim), cbar
xlabelf('units: 1 / %s', sg.units)
im(2, hsino, 'sino'), cbar

% low-dose noise parameters
bi = 1e6; % noise free
k = 0.1924; % low-dose 17 mAs
sig = sqrt(10); % standard variance of electronic noise, a characteristic of CT scanner

% perform low-dose simulation on the projection data
[lsino,yi] = pct_ldsino(bi, hsino, k, sig);

xl = @(x) xlabelf('RMSE = %.3f / %s', rms(col(x - xtrue)), sg.units);

% reconstruct image
tmp = fbp2(sg, ig);
fbp = fbp2(lsino, tmp);
fbpw = fbp2(lsino, tmp, 'window', 'hanning,0.6');

if ~isvar('kappa'), printm 'kappa: try to make resolution approximately uniform'
    wi = yi; % will give 0 weight to any ray where yi=0!
    kappa = sqrt( div0(A' * wi, A' * ones(size(wi))) );
    im(4, kappa), cbar
%     prompt
end


% use local psf to help select beta
if ~isvar('R'), printm 'R'
    f.l2b = 8; % maybe a bit too big, but ok for now
    f.delta = 0.01;
    if f.do_sparse
        f.l2b = 4; % sparse angle
    end
    %	f.pot_arg = {'lange3', f.delta}; % todo: why not as sharp as hyper3?
    f.pot_arg = {'hyper3', f.delta};
    R = Reg1(kappa, 'beta', 2^f.l2b, 'pot_arg', f.pot_arg);
    %	qpwls_psf(A, R, 1, ig.mask, Gdiag(wi), 'loop', 1); % use this to choose beta
end


if ~isvar('xpwls'), printm 'iterative reconstruction'
    if has_mex_jf
        f.niter = 20;
        Ab = Gblock(A, 41); % 41 subsets
        xpwls = pwls_sqs_os(fbp(ig.mask), Ab, lsino, R, 'wi', wi, 'niter', f.niter);
    else
        f.niter = 90;
        xpwls = pwls_pcg1(fbp(ig.mask), A, Gdiag(wi), lsino(:), R, 'niter', f.niter);
    end
    xpwls = ig.embed(xpwls);
    im(4, xpwls, 'PWLS'), cbar
    xl(xpwls)
    %prompt
end

% show results
im plc 2 4
clim = [0 1];
im(1, xtrue, 'high-dose', clim), cbar;
im(2, fbp, 'FBP low-dose', clim), cbar
im(3, fbpw, 'PW low-dose', clim), cbar
im(4, xpwls, 'PWLS low-dose', clim), cbar
xl(xpwls)
im(5, hsino, 'sino'), cbar
im(6, lsino, 'noisy sino'), cbar
