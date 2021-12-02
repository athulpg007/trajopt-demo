%--------------- Multiple-Stage Launch Vehicle Ascent Problem ------------%
clear all; clc
%-------------------------------------------------------------------------%
%--------------- Provide All Physical Data for Problem -------------------%
%-------------------------------------------------------------------------%
marsRadius          = 3389500;       % Mars Radius, [m]
gravParam           = 4.282837e13;   % Mars Gravitational Parameter [m3 s-2]
initialMass         = 39075;         % Lift-off Mass, [kg]; REF: Polsgrove et al. 2016
marsRotRate         = 7.0879e-5;     % Mars rotation rate, [rad/s]; REF: NASA.GOV 
seaLevelDensity     = 0.020;         % Mars surface atm. density, [kg/m3]; REF: NASA.GOV         
densityScaleHeight  = 11100;         % Mars atm. scale height, [m]; REF: NASA.GOV   
g0                  = 9.80665;       % Standard gravity acc. on Earth, [m/s2]

scales.length       = 1;
scales.speed        = 1;
scales.time         = 1;
scales.acceleration = 1;
scales.mass         = 1;
scales.force        = 1;
scales.area         = 1;
scales.volume       = 1;
scales.density      = 1;
scales.gravparam    = 1;

if 0,
scales.length       = marsRadius;
scales.speed        = sqrt(gravParam/scales.length);
scales.time         = scales.length/scales.speed;
scales.acceleration = scales.speed/scales.time;
scales.mass         = initialMass;
scales.force        = scales.mass*scales.acceleration;
scales.area         = scales.length^2;
scales.volume       = scales.area.*scales.length;
scales.density      = scales.mass/scales.volume;
scales.gravparam    = scales.acceleration*scales.length^2;
end

omega               = marsRotRate*scales.time;
auxdata.omegaMatrix = omega*[0 -1 0;1 0 0;0 0 0];
auxdata.mu          = gravParam/scales.gravparam;
auxdata.cd          = 0.8;                              % Vehicle CD
auxdata.sa          = 4*pi/scales.area;                 % Vehicle Ref. Area 
auxdata.rho0        = seaLevelDensity/scales.density;
auxdata.H           = densityScaleHeight/scales.length;
auxdata.Re          = marsRadius/scales.length;  
auxdata.g0          = g0/scales.acceleration;

lat0                = 3.0*pi/180;                     % launch site latitude, [rad]
x0                  = auxdata.Re*cos(lat0);            % x0, in ECI coordinates, [m]
y0                  = 0;                               % y0, in ECI coordinates, [m]
z0                  = auxdata.Re*sin(lat0);            % z0, in ECI coordinates, [m]
r0                  = [x0 y0 z0];                      % r0 = [x0 y0 z0], [m]
v0                  = r0*auxdata.omegaMatrix.';        % v0 = r * omegaMatrix, [m/s]

btFirst  = 238.0/scales.time;  % Max. burn time of first stage, [s];  REF: Polsgrove et al. 2016
btSecond = 331.0/scales.time;  % Max. burn time of second stage, [s]; REF: Polsgrove et al. 2016  

t0 = 0/scales.time;            % Lift-off time = T+0
t1 = 238.0/scales.time;        % First stage cut off = T+238.0
t2 = 569.0/scales.time;        % SECO (max dur=331s) = T+569.0


mTotFirst    = 23473.0/scales.mass;   % First stage total mass, [kg]; REF: Polsgrove et al. 2016
mPropFirst   = 20278.0/scales.mass;   % First stage prop. mass, [kg]; REF: Polsgrove et al. 2016
mDryFirst    = mTotFirst-mPropFirst;  % First stage dry   mass, [kg]; REF: Polsgrove et al. 2016

mTotSecond   = 11756.0/scales.mass;   % Second stage total mass, [kg]; REF: Polsgrove et al. 2016
mPropSecond  = 9377.0 /scales.mass;   % Second stage prop. mass, [kg]; REF: Polsgrove et al. 2016 
mDrySecond   = mTotSecond-mPropSecond;% Second stage dry   mass, [kg]; REF: Polsgrove et al. 2016

mPayload     = 3846.0/scales.mass;      % Payload mass, [kg]; REF: Polsgrove et al. 2016

thrustFirst  = 300.0E3/scales.force;    % First stage thrust, [N]; REF: Polsgrove et al. 2016
thrustSecond = 100.0E3/scales.force;    % Second stage thrust,[N]; REF: Polsgrove et al. 2016

mdotFirst    = mPropFirst/btFirst;                  % First stage mass flow rate, [kg/s]
ispFirst     = thrustFirst/(auxdata.g0*mdotFirst);  % First stage ISP, [s]

mdotSecond   = mPropSecond/btSecond;                % Second stage mass flow rate, [kg/s]
ispSecond    = thrustSecond/(auxdata.g0*mdotSecond);% Second stage ISP, [s]

af       = 3564.5E3/scales.length;      % terminal SMA corr. to 100 x 250 km orbit, [m]
ef       = 0.0210;                      % terminal ecc corr. to 100 x 250 km orbit, [-]
incf     = 3.0*pi/180;                 % terminal inc, [rad]
Omf      = 269.8*pi/180;                % terminal RAAN,[rad]
omf      = 130.5*pi/180;                % terminal AoP, [rad]
nuguess = 0;                            % terminal TA,  [rad]

cosincf = cos(incf);
cosOmf = cos(Omf);
cosomf = cos(omf);
oe = [af ef incf Omf omf nuguess];
[rout,vout] = launchoe2rv(oe,auxdata.mu); % terminal r, v vectors in ECI XYZ
rout = rout';
vout = vout';

m10 = mPayload+mTotSecond+mTotFirst;     % initial mass in phase 1, [kg]
m1f = m10-mdotFirst*t1;                  % final mass in phase 1, [kg]

m20 = m1f-mDryFirst;                     % initial mass in phase 2, [kg]
m2f = m20-(mdotSecond)*(t2-t1);          % final mass in phase 2, [kg]

auxdata.thrustFirst  = thrustFirst;
auxdata.thrustSecond = thrustSecond;

auxdata.ispFirst     = ispFirst;
auxdata.ispSecond    = ispSecond;

rmin = -2*auxdata.Re;
rmax = -rmin;
vmin = -10000/scales.speed;
vmax = -vmin;


%-------------------------------------------------------------------------%
%---------- Provide Bounds and Guess in Each Phase of Problem ------------%
%-------------------------------------------------------------------------%
iphase = 1;
bounds.phase(iphase).initialtime.lower = [t0]; 
bounds.phase(iphase).initialtime.upper = [t0]; 
bounds.phase(iphase).finaltime.lower = [t1]; 
bounds.phase(iphase).finaltime.upper = [t1]; 
bounds.phase(iphase).initialstate.lower = [r0(1:3),v0(1:3),m10];    
bounds.phase(iphase).initialstate.upper = [r0(1:3),v0(1:3),m10];    
bounds.phase(iphase).state.lower = [rmin*ones(1,3),vmin*ones(1,3),m1f];
bounds.phase(iphase).state.upper = [rmax*ones(1,3),vmax*ones(1,3),m10];
bounds.phase(iphase).finalstate.lower = [rmin*ones(1,3),vmin*ones(1,3),m1f]; 
bounds.phase(iphase).finalstate.upper = [rmax*ones(1,3),vmax*ones(1,3),m1f]; 
bounds.phase(iphase).control.lower = -ones(1,3);
bounds.phase(iphase).control.upper = +ones(1,3);
bounds.phase(iphase).path.lower  = 1;
bounds.phase(iphase).path.upper  = 1;
guess.phase(iphase).time = [t0; t1];
guess.phase(iphase).state(:,1) = [r0(1); r0(1)];
guess.phase(iphase).state(:,2) = [r0(2); r0(2)];
guess.phase(iphase).state(:,3) = [r0(3); r0(3)];
guess.phase(iphase).state(:,4) = [v0(1); v0(1)];
guess.phase(iphase).state(:,5) = [v0(2); v0(2)];
guess.phase(iphase).state(:,6) = [v0(3); v0(3)];
guess.phase(iphase).state(:,7) = [m10; m1f];
guess.phase(iphase).control(:,1) = [0; 0];
guess.phase(iphase).control(:,2) = [1; 1];
guess.phase(iphase).control(:,3) = [0; 0];

iphase = 2;
bounds.phase(iphase).initialtime.lower = [t1]; 
bounds.phase(iphase).initialtime.upper = [t1]; 
bounds.phase(iphase).finaltime.lower = [t1]; 
bounds.phase(iphase).finaltime.upper = [t2]; 
bounds.phase(iphase).initialstate.lower = [rmin*ones(1,3),vmin*ones(1,3),m20];
bounds.phase(iphase).initialstate.upper = [rmax*ones(1,3),vmax*ones(1,3),m20];
bounds.phase(iphase).state.lower = [rmin*ones(1,3),vmin*ones(1,3),m2f];
bounds.phase(iphase).state.upper = [rmax*ones(1,3),vmax*ones(1,3),m20];
bounds.phase(iphase).finalstate.lower = [rmin*ones(1,3),vmin*ones(1,3),m2f]; 
bounds.phase(iphase).finalstate.upper = [rmax*ones(1,3),vmax*ones(1,3),m20]; 
bounds.phase(iphase).control.lower = -ones(1,3);
bounds.phase(iphase).control.upper = +ones(1,3);
bounds.phase(iphase).path.lower  = 1;
bounds.phase(iphase).path.upper  = 1;
guess.phase(iphase).time    = [t1; t2];
guess.phase(iphase).state(:,1) = [rout(1) rout(1)];
guess.phase(iphase).state(:,2) = [rout(2) rout(2)];
guess.phase(iphase).state(:,3) = [rout(3) rout(3)];
guess.phase(iphase).state(:,4) = [vout(1) vout(1)];
guess.phase(iphase).state(:,5) = [vout(2) vout(2)];
guess.phase(iphase).state(:,6) = [vout(3) vout(3)];
guess.phase(iphase).state(:,7) = [m20; m2f];
guess.phase(iphase).control(:,1) = [0; 0];
guess.phase(iphase).control(:,2) = [1; 1];
guess.phase(iphase).control(:,3) = [0; 0];

%-------------------------------------------------------------------------%
%------------- Set up Event Constraints That Link Phases -----------------%
%-------------------------------------------------------------------------%
bounds.eventgroup(1).lower = [zeros(1,6), -mDryFirst, 0];
bounds.eventgroup(1).upper = [zeros(1,6), -mDryFirst, 0];

%-------------------------------------------------------------------------%
%----------- Set up Event Constraints That Define Final Orbit ------------%
%-------------------------------------------------------------------------%
bounds.eventgroup(2).lower = [af, ef, incf, Omf, omf];
bounds.eventgroup(2).upper = [af, ef, incf, Omf, omf];

%-------------------------------------------------------------------------%
%-------------- Provide an Initial Mesh in Each Phase --------------------%
%-------------------------------------------------------------------------%
for i=1:2
  meshphase(i).colpoints = 4*ones(1,10);
  meshphase(i).fraction = 0.1*ones(1,10);
end

%-------------------------------------------------------------------------%
%----------- Assemble All Information into Setup Structure ---------------%
%-------------------------------------------------------------------------%
setup.name = 'Launch-Vehicle-Ascent-Problem';
setup.functions.continuous = @launchContinuous;
setup.functions.endpoint = @launchEndpoint;
setup.mesh.phase = meshphase;
setup.nlp.solver = 'snopt';
setup.bounds = bounds;
setup.guess = guess;
setup.auxdata = auxdata;
setup.derivatives.supplier = 'sparseFD';
setup.derivatives.derivativelevel = 'second';
setup.derivatives.dependencies = 'sparseNaN';
setup.scales.method = 'automatic-bounds';
setup.mesh.method = 'hp1';
setup.mesh.tolerance = 1e-6;
setup.method = 'RPMintegration';

%-------------------------------------------------------------------------%
%---------------------- Solve Problem using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
totaltime = tic;
output = gpops2(setup);
totaltime = toc(totaltime);


