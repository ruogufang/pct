function [ out ] = pct_test( in )
%PCT_TEST A "playground" function for trying new syntax and commands
%
%   Kolbeinn Karlsson
%

[t h w] = size(in);

out = zeros(1,t);

for i = 1:t
    out(i) = norm(squeeze(in(i,:,:)),'fro');
end

end

