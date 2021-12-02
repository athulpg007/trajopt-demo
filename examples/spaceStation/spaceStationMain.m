% ----------- Space Station Attitude Control Example -----------%
% This example is taken verbatim from the following reference:  %
% Betts, J. T., Practical Methods for Optimal Control Using     %
% Nonlinear Programming, SIAM Press, Philadelphia, 2009.        %
% --------------------------------------------------------------%
%close all
clear all
clc

%-------------------------------------%
%            Problem Setup            %
%-------------------------------------%
auxdata.J(1,:)   = [2.80701911616e7 4.822509936e5 -1.71675094448e7];
auxdata.J(2,:)   = [auxdata.J(1,2) 9.5144639344e7 6.02604448e4];
auxdata.J(3,:)   = [auxdata.J(1,3) auxdata.J(2,3) 7.6594401336e7];
auxdata.omega0   = [-9.5380685844896e-6 -1.1363312657036e-3 5.3472801108427e-6];
auxdata.omegaorb = 0.06511*pi/180;
auxdata.Jinv     = inv(auxdata.J);
angmommax        = 10000;

%-------------------------------------%
%          Boundary Conditions        %
%-------------------------------------%
t0 = 0; tf = 1800;
omega0 = auxdata.omega0;
r0     = [2.9963689649816e-3 1.5334477761054e-1 3.8359805613992e-3];
h0     = [5000 5000 5000];

%----------------------------------------------------%
% Lower and Upper Limits on Time, State, and Control %
%----------------------------------------------------%
omega1min = -2e-3; omega1max = -omega1min;
omega2min = omega1min; omega2max = omega1max;
omega3min = omega1min; omega3max = omega1max;
omegamin = [omega1min,omega2min,omega3min];
omegamax = [omega1max,omega2max,omega3max];
r1min     = -1; r1max     = +1;
r2min     = -1; r2max     = +1;
r3min     = -1; r3max     = +1;
rmin      = [r1min, r2min, r3min];
rmax      = [r1max, r2max, r3max];
angmombound = 15000;
h1min     = -angmombound; h1max     = angmombound;
h2min     = -angmombound; h2max     = angmombound;
h3min     = -angmombound; h3max     = angmombound;
hmin      = [h1min, h2min, h3min];
hmax      = [h1max, h2max, h3max];
u1min     = -150; u1max = 150;
u2min     = -150; u2max = 150;
u3min     = -150; u3max = 150;
umin      = [u1min, u2min, u3min];     
umax      = [u1max, u2max, u3max];

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf;
bounds.phase.finaltime.upper = tf;
bounds.phase.initialstate.lower = [omega0, r0, h0];
bounds.phase.initialstate.upper = [omega0, r0, h0];
bounds.phase.state.lower = [omegamin,rmin,hmin];
bounds.phase.state.upper = [omegamax,rmax,hmax];
% bounds.phase.finalstate.lower = [omegamin,rmin,hmin];
% bounds.phase.finalstate.upper = [omegamax,rmax,hmax];
bounds.phase.finalstate.lower = [omegamin,rmin,zeros(1,3)];
bounds.phase.finalstate.upper = [omegamax,rmax,zeros(1,3)];
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 10;

bounds.phase.control.lower = umin;
bounds.phase.control.upper = umax;
bounds.phase.path.lower = 0;
bounds.phase.path.upper = angmommax^2;

%-------------------------------------------------------------------------%
%---------- Set up Event Constraints That Terminal Constraints -----------%
%-------------------------------------------------------------------------%
bounds.eventgroup(1).lower = zeros(1,6);
bounds.eventgroup(1).upper = zeros(1,6);

%----------------------%
% Set up Initial Guess %
%----------------------%
tGuess               = [t0; tf];
omegaGuess           = [omega0; omega0];
rGuess               = [r0; r0];
hGuess               = [h0; h0];
uGuess               = [zeros(1,3); zeros(1,3)];
guess.phase.state    = [omegaGuess, rGuess, hGuess];
guess.phase.control  = [uGuess];
guess.phase.time     = tGuess;
guess.phase.integral = 10;

%---------------------%
% Set up Initial Mesh %
%---------------------%
NN = 1;
meshphase.colpoints = 4*ones(1,NN);
meshphase.fraction = 1/NN*ones(1,NN);

%-------------------------------------%
%          Problem Setup              %
%-------------------------------------%
setup.name = 'Space-Station-Attitude-Control';
setup.functions.continuous = @spaceStationContinuous;
setup.functions.endpoint   = @spaceStationEndpoint;
setup.auxdata = auxdata;
setup.mesh.phase = meshphase;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'ipopt';
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.derivatives.stepsize = 1e-4;
setup.scales.method = 'automatic-bounds';
setup.method = 'RPMintegration';
setup.mesh.method = 'hp1';
setup.mesh.tolerance = 1e-6;
setup.mesh.colpointsmin = 4;
setup.mesh.colpointsmax = 16;

%----------------------------------%
%   Solve Problem Using GPOPS-II   %
%----------------------------------%
output = gpops2(setup);
