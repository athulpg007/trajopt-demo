N = 20;
thetavec = linspace(0,pi,N);
phivec = linspace(0,2*pi,2*N);
[th, ph] = meshgrid(thetavec,phivec);

R = 6378145.0*ones(size(th)); % should be your R(theta,phi) surface in general

x = R.*sin(th).*cos(ph);
y = R.*sin(th).*sin(ph);
z = R.*cos(th);




solution = output.result.solution;

X1=solution.phase(1).state(:,1);
Y1=solution.phase(1).state(:,2);
Z1=solution.phase(1).state(:,3);

X2=solution.phase(2).state(:,1);
Y2=solution.phase(2).state(:,2);
Z2=solution.phase(2).state(:,3);

X3=solution.phase(3).state(:,1);
Y3=solution.phase(3).state(:,2);
Z3=solution.phase(3).state(:,3);

X4=solution.phase(4).state(:,1);
Y4=solution.phase(4).state(:,2);
Z4=solution.phase(4).state(:,3);

plot3(X1,Y1,Z1)
plot3(X2,Y2,Z2)
plot3(X3,Y3,Z3)
plot3(X4,Y4,Z4)

hold('on')

surf(x,y,z);
axis vis3d