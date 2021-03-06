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
t{3} = solution.phase(3).time;
rad3 = solution.phase(3).state(:,1:3);
alt{3} = (sqrt(dot(rad3,rad3,2))-auxdata.Re)*scales.length/1000;
velocity3 = solution.phase(3).state(:,4:6);
speed{3} = sqrt(dot(velocity3,velocity3,2));
mass{3} = solution.phase(3).state(:,7);
control{3} = solution.phase(3).control;
t{4} = solution.phase(4).time;
rad4 = solution.phase(4).state(:,1:3);
alt{4} = (sqrt(dot(rad4,rad4,2))-auxdata.Re)*scales.length/1000;
speed{4} = solution.phase(4).state(:,4:6);
velocity4 = solution.phase(4).state(:,4:6);
speed{4} = sqrt(dot(velocity4,velocity4,2));
mass{4} = solution.phase(4).state(:,7);
control{4} = solution.phase(4).control;

% Plot the Altitude in Each Phase
figure(1)
pp = plot(t{1},alt{1},'-o',t{2},alt{2},'-o',t{3},alt{3},'-o',t{4},alt{4},'-o');
xl = xlabel('time (s)');
yl = ylabel('altitude (km)');
ll = legend('Phase 1','Phase 2','Phase 3','Phase 4','Location','SouthEast');
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
pp = plot(t{1},speed{1},'-o',t{2},speed{2},'-o',t{3},speed{3},'-o',t{4},speed{4},'-o');
xl = xlabel('time (s)');
yl = ylabel('speed (m/s)');
ll = legend('Phase 1','Phase 2','Phase 3','Phase 4','Location','SouthEast');
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
pp = plot(t{1},control{1},'-o',t{2},control{2},'-o',t{3},control{3},'-o',t{4},control{4},'-o');
xl = xlabel('time (s)');
yl = ylabel('control');
ll = legend('Phase 1','Phase 2','Phase 3','Phase 4','Location','SouthWest');
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
pp = plot(t{1},mass{1},'-o',t{2},mass{2},'-o',t{3},mass{3},'-o',t{4},mass{4},'-o');
xl = xlabel('time (s)');
yl = ylabel('mass');
ll = legend('Phase 1','Phase 2','Phase 3','Phase 4','Location','NorthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.5);
grid on
print -depsc2 launchMass.eps
print -dpng launchMass.png
