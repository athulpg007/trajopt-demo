iconvert = 0;
solution = output.result.solution;
sos = load('sosSolution.dat');
sosTime = sos(:,1);
sosState = sos(:,2:10);
sosControl = sos(:,11:13);
time  = solution.phase.time;
state = solution.phase.state;
control = solution.phase.control;

ifig = 1;
figure(ifig);
pp = plot(sosTime,sosState(:,1)*1e4,'-bo',time,state(:,1)*1e4,'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\omega_1(t)\times 10^{-4}$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('omega1 e-4');
ll = legend('SOS','GPOPS-II','Location','SouthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState1.eps
if iconvert,
  eval('!eps2png spaceStationState1.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosState(:,2)*1e4,'-bo',time,state(:,2)*1e4,'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\omega_1(t)\times 10^{-4}$','Interpreter','LaTeX');
end;
xl = xlabel('time');
yl = ylabel('omega2 e-4');
ll = legend('SOS','GPOPS-II','Location','Best');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState2.eps
if iconvert,
  eval('!eps2png spaceStationState2.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosState(:,3)*1e4,'-bo',time,state(:,3)*1e4,'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\omega_3(t)\times 10^{-4}$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('omega3 e-4');
ll = legend('SOS','GPOPS-II','Location','NorthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState3.eps
if iconvert,
  eval('!eps2png spaceStationState3.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosState(:,4),'-bo',time,state(:,4),'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$r_1(t)$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('r1');
ll = legend('SOS','GPOPS-II','Location','NorthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState4.eps
if iconvert,
  eval('!eps2png spaceStationState4.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosState(:,5),'-bo',time,state(:,5),'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$r_2(t)$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('r2');
ll = legend('SOS','GPOPS-II','Location','NorthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState5.eps
if iconvert,
  eval('!eps2png spaceStationState5.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosState(:,6),'-bo',time,state(:,6),'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$r_3(t)$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('r3');
ll = legend('SOS','GPOPS-II','Location','NorthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState6.eps
if iconvert,
  eval('!eps2png spaceStationState6.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosState(:,7)/1000,'-bo',time,state(:,7)/1000,'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$h_1(t)\times 1000$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('h1');
ll = legend('SOS','GPOPS-II','Location','NorthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState7.eps
if iconvert,
  eval('!eps2png spaceStationState7.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosState(:,8)/1000,'-bo',time,state(:,8)/1000,'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$h_2(t)\times 1000$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('h2');
ll = legend('SOS','GPOPS-II','Location','NorthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState8.eps
if iconvert,
  eval('!eps2png spaceStationState8.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosState(:,9)/1000,'-bo',time,state(:,9)/1000,'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$h_3(t)\times 1000$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('h3');
ll = legend('SOS','GPOPS-II','Location','SouthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationState9.eps
if iconvert,
  eval('!eps2png spaceStationState9.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosControl(:,1),'-bo',time,control(:,1),'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_1(t)$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('u1');
ll = legend('SOS','GPOPS-II','Location','SouthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationControl1.eps
if iconvert,
  eval('!eps2png spaceStationControl1.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosControl(:,2),'-bo',time,control(:,2),'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_2(t)$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('u2');
ll = legend('SOS','GPOPS-II','Location','SouthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationControl2.eps
if iconvert,
  eval('!eps2png spaceStationControl2.eps');
end;

ifig = ifig+1;
figure(ifig);
pp = plot(sosTime,sosControl(:,3),'-bo',time,control(:,3),'-rd');
set(pp,'LineWidth',1.25,'MarkerSize',6);
if 0,
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_3(t)$','Interpreter','LaTeX');
end
xl = xlabel('time');
yl = ylabel('u3');
ll = legend('SOS','GPOPS-II','Location','SouthEast');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',20);
set(gca,'FontSize',18,'FontName','Times');
grid on;
print -depsc2 spaceStationControl3.eps
if iconvert,
  eval('!eps2png spaceStationControl3.eps');
end;

