
% Generate Ca and R. Compute c from Ca and R.
% Set \Delta t and F_t to 1 for simplicity.
% Bolus delay d_t = 1 s

% Ground truth Ca
Ca = [1 0 0 0;
    2 1 0 0;
    0 2 1 0
    0 0 2 1];

% delayed Ca
Cad = [0 0 0 0;
	1 0 0 0;
    2 1 0 0
    0 2 1 0];

R = [1, 0.5 0 0]';	% This is our ground-truth R

% Normalize R to 1
% R = R / sum(R);

% Compute c from Ca and R
c = Ca * R;

% Estimate R by regular SVD
[U, S, V] = svd(Cad);
r_est1 = V * pinv(S) * U' * c;

% Zero-pad Ca and c
Cap = [0 0 0 0 0 0 2 1;
    1 0 0 0 0 0 0 2;
    2 1 0 0 0 0 0 0;
    0 2 1 0 0 0 0 0;
    0 0 2 1 0 0 0 0;
    0 0 0 2 1 0 0 0
    0 0 0 0 2 1 0 0
    0 0 0 0 0 2 1 0];

cp = [c; zeros(length(c), 1)];

% Re-estimate R with 
[Up, Sp, Vp] = svd(Cap);
r_est2 = Vp * pinv(Sp) * Up' * cp;

% Print results to screen
disp('Ground truth R:')
disp(R')
disp('SVD-estimated R:')
disp(r_est1')
disp('Block-circulant SVD-estimated R:')
disp(r_est2')
disp('As you can see, R(N+1:L) is non-zero.')
disp('But I think one should discard the R(N+1:L) anyway.')


% In this case, the block-circulant SVD gives 
% a worse result, since c is a noise-free 
% convolution of 



	
