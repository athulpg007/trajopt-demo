%-------------------------------------------------------------------------%
%                             Extract Solution                            %
%-------------------------------------------------------------------------%
SOCS = load('hyperSensitiveSolutionSOCS10K.dat');
solution = output.result.solution;
time = solution.phase(1).time;
state = solution.phase(1).state;
control = solution.phase(1).control;
for i=1:length(output.meshhistory);
  mesh(i).points = [0 cumsum(output.meshhistory(i).result.setup.mesh.phase.fraction)];
  mesh(i).iteration = i*ones(size(mesh(i).points));
end;

%-------------------------------------------------------------------------%
%                              Plot Solution                              %
%-------------------------------------------------------------------------%
figure(1);
pp = plot(SOCS(:,1),SOCS(:,2),'-bo',time,state,'-rd');
set(pp,'MarkerSize',6,'LineWidth',1.25);
xl = xlabel('time');
yl = ylabel('state');
ll = legend('SOCS','GPOPS-II');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(ll,'FontSize',20);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 hyperSensitiveState.eps
print -dpng hyperSensitiveState.png

figure(2);
pp = plot(SOCS(:,1),SOCS(:,3),'-bo',time,control,'-rd');
set(pp,'MarkerSize',6,'LineWidth',1.25);
xl = xlabel('time');
yl = ylabel('control');
ll = legend('SOCS','GPOPS-II');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(ll,'FontSize',20);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 hyperSensitiveControl.eps
print -dpng hyperSensitiveControl.png

figure(3);
pp = plot(SOCS(:,1),SOCS(:,2),'-bo',time,state,'-rd');
set(pp,'MarkerSize',6,'LineWidth',1.25);
xl = xlabel('time');
yl = ylabel('state');
ll = legend('SOCS','GPOPS-II');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(ll,'FontSize',20);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
axis([0,25,-0.2,1.6]);
print -depsc2 hyperSensitiveStateInitialLayer.eps
print -dpng hyperSensitiveStateInitialLayer.png

figure(4);
pp = plot(SOCS(:,1),SOCS(:,2),'-bo',time,state,'-rd');
set(pp,'MarkerSize',6,'LineWidth',1.25);
xl = xlabel('time');
yl = ylabel('state');
ll = legend('SOCS','GPOPS-II');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(ll,'FontSize',20);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
axis([9975,10000,-0.2,1.6]);
print -depsc2 hyperSensitiveStateFinalLayer.eps
print -dpng hyperSensitiveStateFinalLayer.png

figure(5);
pp = plot(SOCS(:,1),SOCS(:,3),'-bo',time,control,'-rd');
set(pp,'MarkerSize',6,'LineWidth',1.25);
xl = xlabel('time');
yl = ylabel('control');
ll = legend('SOCS','GPOPS-II');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(ll,'FontSize',20);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
axis([0,25,-1,1]);
print -depsc2 hyperSensitiveControlInitialLayer.eps
print -dpng hyperSensitiveControlInitialLayer.png

figure(6);
pp = plot(SOCS(:,1),SOCS(:,3),'-bo',time,control,'-rd');
set(pp,'MarkerSize',6,'LineWidth',1.25);
xl = xlabel('time');
yl = ylabel('control');
ll = legend('SOCS','GPOPS-II');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(ll,'FontSize',20);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
axis([9975,10000,-1,8]);
print -depsc2 hyperSensitiveControlFinalLayer.eps
print -dpng hyperSensitiveControlFinalLayer.png

figure(7);
for i=1:length(mesh);
  pp = plot(mesh(i).points*tf,mesh(i).iteration,'bo');
  set(pp,'LineWidth',1.25);
  hold on;
end;
xl = xlabel('Mesh Point Location (Fraction of Interval)');
yl = ylabel('Mesh Iteration');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'YTick',0:length(mesh),'FontSize',16);
grid on;
print -depsc2 hyperSensitiveMeshRefinement.eps
print -dpng hyperSensitiveMeshRefinement.png
