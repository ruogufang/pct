function [ out ] = thfunc(dt,Ca,F,E,Ve,Tc)
%THFUNC Tissue Homogeneity function 
%
%   USAGE: OUT = THFUNC
%   More to come
%
%   Is used in the optimization in pct_acth and pct_acthp.

%Preliminaries
L = length(Ca); 

k = E*F/Ve;
delay = Tc/dt;
t = 0:dt:((L-delay)*dt-dt);

%Intravascular component of R
Ri = [ones(delay,1); zeros(L-delay,1)];

%Extravascular component of R
if Ve ~= 0
    Re = [zeros(delay,1); E*exp(-k*t)'];
else
    Re = zeros(size(Ri));
end
%R
R = Ri + Re;

Q = conv(F*R,Ca);

out = Q(1:L);

if max(abs(out)) > 1e20
    out(:) = 0;
    


end