function [Ca,Cc,Rc] = pct_circ(AIF,C,R,m)
%PCT_CIRC creates block-circulant version of Ca, C and R
%
%   Ruogu Fang Revised 09/12/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  [Ca,Cc,Rc] = PCT_CIRC(AIF,C,R,M);
%
%   PRE:
%       AIF     - Arterial input function [T x 1]
%       C       - CTP time attenuation curves (TAC) [T x Y x X]
%       R       - Input residue functions [T x Y x X]
%       m       - block-circulant factor (extend time T to m*T)
%
%   POST:
%       Ca      - Circulant matrix of AIF [mT x mT]
%       Cc      - Extended CAT [mT x Y x X]
%       Rc      - Extended IRF [mT x Y x X]
%
% PCT_CIRC creates block-circulant version of arterial input function
% matrix Ca, and extends C and R by zero padding on the time axis.
%

if m == 1 % no block-circulant
    Ca = tril(toeplitz(AIF));
    Cc = C;
    Rc = R;
elseif m >= 2 % block-circulant
    T = size(C,1);
    if size(AIF,1)<size(AIF,2)
        AIF = AIF';
    end
    AIF = padarray(AIF,(m-1)*T,'post');
    Cc = padarray(C,(m-1)*T,'post');
    Rc = padarray(R,(m-1)*T,'post');
    
    %Create the block-circulant matrix
    Ca = sparse(gallery('circul',AIF)');
end
end