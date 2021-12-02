% Modified from GPOPS II launch example code.
% Edited by: Athul P. G, apradee@purdue.edu, 5/6/2018
% Puropse: Define problem state equations for each phase.

function phaseout = launchContinuous(input)

%---------------------%
% Dynamics in Phase 1 / First Stage 1 %
%---------------------%
t1 = input.phase(1).time;    % t1 = phase 1 time 
x1 = input.phase(1).state;   % x1 = phase 1 state
u1 = input.phase(1).control; % u1 = phase 1 control

% extract state variables r = X,Y,Z; v = VX,VY,VZ, m
% for phase 1
r1 = x1(:,1:3);             % r1 = phase 1 XYZ
v1 = x1(:,4:6);             % v1 = phase 1 VXVYVZ
m1 = x1(:,7);               % m1 = phase 1 m

% compute radial distance r from sqrt(x^2+y^2+z^2)
rad1 = sqrt(sum(r1.*r1,2));
% compute relative vel vrel = v - \omega x r
omegaMatrix = input.auxdata.omegaMatrix;
omegacrossr = r1*omegaMatrix.';
vrel1 = v1-omegacrossr;
% compute speed rel = |vrel|
speedrel1 = sqrt(sum(vrel1.*vrel1,2));
% compute altitude h
h1 = rad1-input.auxdata.Re;
% compute density at altitude h
rho1 = input.auxdata.rho0*exp(-h1/input.auxdata.H);

bc1  = (rho1./(2*m1)).*input.auxdata.sa*input.auxdata.cd;
bcspeed1 = bc1.*speedrel1;
bcspeedmat1 = repmat(bcspeed1,1,3);
Drag1 = -bcspeedmat1.*vrel1;
muoverradcubed1 = input.auxdata.mu./rad1.^3;

muoverradcubedmat1 = [muoverradcubed1 muoverradcubed1 muoverradcubed1];
grav1 = -muoverradcubedmat1.*r1;

TFirst1 = input.auxdata.thrustFirst*ones(size(t1));
TTot1   = TFirst1;

m2dot1  = -TFirst1./(input.auxdata.g0*input.auxdata.ispFirst);
mdot1   = m2dot1;

path1 = [sum(u1.*u1,2)];
Toverm1 = TTot1./m1;
Tovermmat1 = [Toverm1 Toverm1 Toverm1];
thrust1 = Tovermmat1.*u1;
rdot1 = v1;
vdot1 = thrust1+Drag1+grav1;
phaseout(1).dynamics = [rdot1 vdot1 mdot1];
phaseout(1).path = path1;

%---------------------%
% Dynamics in Phase 2 / Second Stage %
%---------------------%
t2 = input.phase(2).time;
x2 = input.phase(2).state;
u2 = input.phase(2).control;
r2 = x2(:,1:3);
v2 = x2(:,4:6);
m2 = x2(:,7);

rad2 = sqrt(sum(r2.*r2,2));
omegaMatrix = input.auxdata.omegaMatrix;
omegacrossr = r2*omegaMatrix.';
vrel2 = v2-omegacrossr;
speedrel2 = sqrt(sum(vrel2.*vrel2,2));
h2 = rad2-input.auxdata.Re;
rho2 = input.auxdata.rho0*exp(-h2/input.auxdata.H);
bc2  = (rho2./(2*m2)).*input.auxdata.sa*input.auxdata.cd;
bcspeed2 = bc2.*speedrel2;
bcspeedmat2 = repmat(bcspeed2,1,3);
Drag2 = -bcspeedmat2.*vrel2;
muoverradcubed2 = input.auxdata.mu./rad2.^3;
muoverradcubedmat2 = [muoverradcubed2 muoverradcubed2 muoverradcubed2];
grav2 = -muoverradcubedmat2.*r2;


TSecond2 = input.auxdata.thrustSecond*ones(size(t2)); 
TTot2    = TSecond2;

m2dot2 = -TSecond2./(input.auxdata.g0*input.auxdata.ispSecond);
mdot2  = m2dot2;  


path2 = [sum(u2.*u2,2)];
Toverm2 = TTot2./m2;
Tovermmat2 = [Toverm2 Toverm2 Toverm2];
thrust2 = Tovermmat2.*u2;
rdot2 = v2;
vdot2 = thrust2+Drag2+grav2;
phaseout(2).dynamics = [rdot2 vdot2 mdot2];
phaseout(2).path = path2;

