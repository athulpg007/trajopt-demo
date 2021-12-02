% Modified from GPOPS II launch example code.
% Edited by: Athul P. G, apradee@purdue.edu, 5/6/2018
% Puropse: Coordinate transformation from orbital elements to position
% and velocity vectors in ECI XYZ frame.

function [ri,vi] = launchoe2rv(oe,mu)
%{
Transform Keplerian orbital elements a,e,i,O,o,nu
into to R = X,Y,Z ; V = VX, VY, VZin ECI XYZ frame.

Ref: Bate, 1971 / Astrodynamics 532 Qualifying Exams Questions.

Inputs: oe (array of length 6), mu
Output: r, v (array of length 3)
%}
% extract a, e, i, O, o, nu from input array oe
a=oe(1); e=oe(2); i=oe(3); Om=oe(4); om=oe(5); nu=oe(6);
p = a*(1-e*e);
r = p/(1+e*cos(nu));
% compute r and v vector in EPH frame (perifocal)
rv = [r*cos(nu); r*sin(nu); 0];
vv = sqrt(mu/p)*[-sin(nu); e+cos(nu); 0];
% compute elements of transformation matrix  
cO = cos(Om);  sO = sin(Om);
co = cos(om);  so = sin(om);
ci = cos(i);   si = sin(i);
% transformation matrix R from EPH to ECI XYZ as defined in Bate
R  = [cO*co-sO*so*ci  -cO*so-sO*co*ci  sO*si;
      sO*co+cO*so*ci  -sO*so+cO*co*ci -cO*si;
      so*si            co*si           ci];
% transform rv, vv to ECI XYZ coordinates, return
ri = R*rv;
vi = R*vv;
