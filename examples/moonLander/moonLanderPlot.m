%-------------------------------------------------------------------------%
%                             Extract Solution                            %
%-------------------------------------------------------------------------%
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
pp = plot(time,state(:,1),'-o');
xl = xlabel('time');
yl = ylabel('Altitude');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 moonLanderAltitude.eps
print -dpng moonLanderAltitude.png

figure(2);
pp = plot(time,state(:,2),'-o');
xl = xlabel('time');
yl = ylabel('Velocity');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 moonLanderVelocity.eps
print -dpng moonLanderVelocity.png

figure(3);
pp = plot(time,control,'-o');
xl = xlabel('time');
yl = ylabel('control');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
grid on
print -depsc2 moonLanderControl.eps
print -dpng moonLanderControl.png

figure(4);
tf = time(end);
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
print -depsc2 moonLanderMeshRefinement.eps
print -dpng moonLanderMeshRefinement.png
