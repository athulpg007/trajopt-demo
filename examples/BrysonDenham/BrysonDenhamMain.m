%------------------------- Bryson-Denham Problem -------------------------%
% This problem is taken from the following reference:                     %
% Bryson, A. E. and Ho, Y-C, "Applied Optimal Control:  Optimization,     %
% Estimation, and Control," Hemisphere Publishing, 1975.                  %
%-------------------------------------------------------------------------%
clear all
close all
clc

%-------------------------------------------------------------------------%
%------------- Provide and Set Up All Bounds for Problem -----------------%
%-------------------------------------------------------------------------%
L = 1/9;
t0                              = 0;     
tf                              = 1;
x0                              = 0;
xf                              = 0;
xmin                            = 0;
xmax                            = L;
v0                              = 1;
vf                              = -v0;
vmin                            = -10;
vmax                            = 10;
umin                            = -10;
umax                            = 5;
bounds.phase.initialtime.lower  = t0;
bounds.phase.initialtime.upper  = t0;
bounds.phase.finaltime.lower    = tf;
bounds.phase.finaltime.upper    = tf;
bounds.phase.initialstate.lower = [x0,v0];
bounds.phase.initialstate.upper = [x0,v0];
bounds.phase.state.lower        = [xmin,vmin];
bounds.phase.state.upper        = [xmax,vmax];
bounds.phase.finalstate.lower   = [xf,vf];
bounds.phase.finalstate.upper   = [xf,vf];
bounds.phase.control.lower      = [umin];
bounds.phase.control.upper      = [umax];
bounds.phase.integral.lower     = 0;
bounds.phase.integral.upper     = 10;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
tGuess               = [t0;tf];
xGuess               = [x0;xf];
vGuess               = [v0;vf];
uGuess               = [umin;0];
guess.phase.time     = tGuess;
guess.phase.state    = [xGuess,vGuess];
guess.phase.control  = uGuess;
guess.phase.integral = 0;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp1';
mesh.tolerance       = 1e-8;
mesh.maxiteration    = 45;
mesh.colpointsmax    = 4;
mesh.colpointsmin    = 10;
N                    = 10;
mesh.phase.colpoints = 4*ones(1,N);
mesh.phase.fraction  =  ones(1,10)/N;


%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                        = 'Bryson-Denham-State-Constrained-Double-Integrator';
setup.functions.continuous        = @BrysonDenhamContinuous;
setup.functions.endpoint          = @BrysonDenhamEvents;
setup.bounds                      = bounds;
setup.guess                       = guess;
setup.nlp.solver                  = 'ipopt';
setup.derivatives.supplier        = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.method                      = 'RPMintegration';
setup.mesh                        = mesh;

%-------------------------------------------------------------------------%
%---------------------- Solve Problem Using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);

