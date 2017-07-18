function [ out ] = pct_objfunc(tac,dt,aif,x1,x2,x3,tcc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

out = sqrt( sum( (tac - thfunc(dt,aif,x1,x2,x3,tcc)).^2) / length(aif));

if out == Inf || out == -Inf || isnan(out)
    out = 0;
end

end

