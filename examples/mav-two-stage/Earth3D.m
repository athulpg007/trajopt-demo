N = 20;
thetavec = linspace(0,pi,N);
phivec = linspace(0,2*pi,2*N);
[th, ph] = meshgrid(thetavec,phivec);

R = 3389500*ones(size(th)); % should be your R(theta,phi) surface in general

x = R.*sin(th).*cos(ph);
y = R.*sin(th).*sin(ph);
z = R.*cos(th);


axis('equal')

solution = output.result.solution;

X1=solution.phase(1).state(:,1);
Y1=solution.phase(1).state(:,2);
Z1=solution.phase(1).state(:,3);

X2=solution.phase(2).state(:,1);
Y2=solution.phase(2).state(:,2);
Z2=solution.phase(2).state(:,3);

hold('on')

plot3(X1,Y1,Z1,'yo')
plot3(X2,Y2,Z2,'bo')

xorb = zeros(121);
yorb = zeros(121);
zorb = zeros(121);

for nu = 1:1:121
    XYZ = launchoe2rv([af,ef,incf,Omf,omf,(nu/101)*2*pi],gravParam);
    xorb(nu) = XYZ(1,1);
    yorb(nu) = XYZ(2,1);
    zorb(nu) = XYZ(3,1);
end



plot3(xorb,yorb,zorb,'g-','LineWidth',3.0)

planeX = [1 -1 -1  1]*5E6;
planeY = [1  1 -1 -1]*5E6;
patch(planeX,planeY,'red','FaceAlpha',0.2)

quiver3(0.0,0.0,0.0,5.0E6,0,0,'LineWidth',2.0)
quiver3(0.0,0.0,0.0,0.0E6,5.0E6,0,'LineWidth',2.0)
quiver3(0.0,0.0,0.0,0.0E6,0,5.0E6,'LineWidth',2.0)


hSurface=surf(x,y,z);
set(hSurface,'FaceColor',[0.9100 0.4100 0.1700],'FaceAlpha',1.0);
axis vis3d