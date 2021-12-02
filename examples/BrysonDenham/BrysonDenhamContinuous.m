function phaseout = BrysonDenhamContinuous(input)

v = input.phase.state(:,2);
u = input.phase.control;

dx                 = v;
dv                 = u;
integrand          = (u.^2)./2;
phaseout.dynamics  = [dx,dv];
phaseout.integrand = integrand;
