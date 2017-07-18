function B = pct_mtxdiag(A)
% B = PCT_MTXDIAG(A) diagonizes the columns of matrix A 
%
% e.g. A = [1 2 3
%           4 5 6]
%
%      B = [1 0 0
%           4 0 0
%           0 2 0
%           0 5 0
%           0 0 3
%           0 0 6]
%
%   Ruogu Fang 09/12/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  B = PCT_MTXDIAG(A);
%
%   PRE:
%       A       - Input matrix [M x N]
%
%   POST:
%       B       - Output sparse matrix [MN x N]
%
% B = PCT_MTXDIAG(A) diagonizes the columns of matrix A [M x N] to matrix B [MN x N] by
% putting the columns of A in the diagnolal position of the new matrix

[M,N] = size(A);
idx_y = (1:M*N)';
idx_x = reshape(repmat(1:N,M,1),[],1);
s = reshape(A,[],1);
B = sparse(idx_y,idx_x,s);

end