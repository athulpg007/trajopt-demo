% Modified from GPOPS II launch example code.
% Edited by: Athul P. G, apradee@purdue.edu, 5/6/2018
% Puropse: Coordinate transformation from position
% and velocity vectors in ECI XYZ frame to orbital elements.

function oe = launchrv2oe(rv,vv,mu);

%{
Transform rv, vv into Keplerian elements.

Ref: Benson, Ph.D Thesis, MIT

Inputs: oe (array of length 6), mu
Output: r, v (array of length 3)
%}
% K = unit vector in ECI Z direction
K  = [0;0;1];
% hv = r x v
hv = cross(rv,vv);
% n = K x h
nv = cross(K,hv);
% n = |n|
n  = sqrt(nv.'*nv);
%h2 = |h|^2
h2 = (hv.'*hv);
%v2 = |v|^2
v2 = (vv.'*vv);
%r = |r|
r  = sqrt(rv.'*rv);
% ev = eccentricity vector 
ev = 1/mu *( (v2-mu/r)*rv - (rv.'*vv)*vv );
% p = h^2 / mu
p  = h2/mu;
% e = |ev|
e  = sqrt(ev.'*ev);
% a = SMA
a  = p/(1-e*e);
% i = inclination
i  = acos(hv(3)/sqrt(h2));
% Om = RAAN
Om = acos(nv(1)/n);
if (nv(2)<0-eps),
  Om = 2*pi-Om;
end;
% om = AoP
om = acos(nv.'*ev/n/e);
if (ev(3)<0),
  om = 2*pi-om;
end;
nu = acos(ev.'*rv/e/r);
if (rv.'*vv<0),
  nu = 2*pi-nu;
end;
% return orbital elememets array oe
oe = [a; e; i; Om; om; nu];
