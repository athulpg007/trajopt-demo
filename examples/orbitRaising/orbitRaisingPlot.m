%------------------------------%
% Extract Solution from Output %
%------------------------------%
solution = output.result.solution;
time = solution.phase(1).time;
r         = solution.phase(1).state(:,1);
theta     = solution.phase(1).state(:,2);
vr        = solution.phase(1).state(:,3);
vtheta    = solution.phase(1).state(:,4);
u1        = solution.phase(1).control(:,1);
u2        = solution.phase(1).control(:,2);

%---------------%
% Plot Solution %
%---------------%
figure(1)
pp = plot(time,r,'-o');
xl = xlabel('time (s)');
yl = ylabel('radius');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 orbitRaisingRadius.eps
print -dpng orbitRaisingRadius.png

figure(2)
pp = plot(time,theta,'-o');
xl = xlabel('time (s)');
yl = ylabel('theta');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 orbitRaisingTheta.eps
print -dpng orbitRaisingTheta.png

figure(3)
pp = plot(time,vr,'-o');
xl = xlabel('time (s)');
yl = ylabel('vr');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 orbitRaisingVr.eps
print -dpng orbitRaisingVr.png

figure(4)
pp = plot(time,vtheta,'-o');
xl = xlabel('time (s)');
yl = ylabel('vtheta');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25,'MarkerSize',8);
grid on
print -depsc2 orbitRaisingVtheta.eps
print -dpng orbitRaisingVtheta.png

figure(4)
pp = plot(time,u1,'-o',time,u2,'-o');
set(pp,'LineWidth',1.25,'MarkerSize',8);
xl = xlabel('time (s)');
yl = ylabel('control');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
grid on
print -depsc2 orbitRaisingVtheta.eps
print -dpng orbitRaisingVtheta.png
