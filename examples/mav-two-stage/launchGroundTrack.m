solution = output.result.solution;

X1=solution.phase(1).state(:,1);
Y1=solution.phase(1).state(:,2);
Z1=solution.phase(1).state(:,3);

X2=solution.phase(2).state(:,1);
Y2=solution.phase(2).state(:,2);
Z2=solution.phase(2).state(:,3);


lon1 = atan2(Y1,X1)*180/pi - 80.6077;
lon2 = atan2(Y2,X2)*180/pi - 80.6077;

lat1 = atan2(Z1, sqrt(X1.^2+Y1.^2))*180/pi;
lat2 = atan2(Z2, sqrt(X2.^2+Y2.^2))*180/pi;

lat = cat(1, lat1, lat2);
lon = cat(1, lon1, lon2);

wmline(lat,lon)



