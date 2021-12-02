% Extract the altitude, speed, and control in each phase of the problem. 
solution = output.result.solution;
t{1} = solution.phase(1).time;
rad1 = solution.phase(1).state(:,1:3);
alt{1} = (sqrt(dot(rad1,rad1,2))-auxdata.Re)*scales.length/1000;
velocity1 = solution.phase(1).state(:,4:6);
speed{1} = sqrt(dot(velocity1,velocity1,2));
mass{1} = solution.phase(1).state(:,7);
control{1} = solution.phase(1).control;
t{2} = solution.phase(2).time;
rad2 = solution.phase(2).state(:,1:3);
alt{2} = (sqrt(dot(rad2,rad2,2))-auxdata.Re)*scales.length/1000;
velocity2 = solution.phase(2).state(:,4:6);
speed{2} = sqrt(dot(velocity2,velocity2,2));
mass{2} = solution.phase(2).state(:,7);
control{2} = solution.phase(2).control;

% Plot the Altitude in Each Phase
figure(1)
pp = plot(t{1},alt{1},'-o',t{2},alt{2},'-o');
xl = xlabel('time (s)');
yl = ylabel('altitude (km)');
ll = legend('Phase 1','Phase 2','Location','SouthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 launchAltitude.eps
print -dpng launchAltitude.png

% Plot the Speed in Each Phase
figure(2)
pp = plot(t{1},speed{1},'-o',t{2},speed{2},'-o');
xl = xlabel('time (s)');
yl = ylabel('speed (m/s)');
ll = legend('Phase 1','Phase 2','Location','SouthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.5);
grid on
print -depsc2 launchSpeed.eps
print -dpng launchSpeed.png

% Plot the Control in Each Phase
figure(3)
pp = plot(t{1},control{1},'-o',t{2},control{2},'-o');
xl = xlabel('time (s)');
yl = ylabel('control');
ll = legend('Phase 1','Phase 2','Location','SouthWest');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.5);
grid on
print -depsc2 launchControl.eps
print -dpng launchControl.png

% Plot the mass in Each Phase
figure(4)
pp = plot(t{1},mass{1},'-o',t{2},mass{2},'-o');
xl = xlabel('time (s)');
yl = ylabel('mass');
ll = legend('Phase 1','Phase 2','Location','NorthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.5);
grid on
print -depsc2 launchMass.eps
print -dpng launchMass.png
