clc;
close all;
clear all;

xi = -10;
xf = 10;
dx = 0.1;
er = 0.2;

a2 = 0.5;
a1 = 1.57;
a0 = -2;

x = xi:dx:xf;
y = a2*x.^2+a1*x+a0;
rg = y(end)-y(1);

nm = 25;
px = sort((xf - xi)*rand(nm,1)+xi);
py = sort(a2*px.^2+a1*px+a0+2*rg*er*rand(nm,1)-rg*er);

x2 = px.*px;
x3 = x2.*px;

Sx = ones(1,nm)*px;
Sy = ones(1,nm)*py;
Sxy = px'*py;
Sx2 = px'*px;
Sx3 = x2'*px;
Sx4 = x3'*px;
Sx2y = x2'*py;

aux = (Sx4*(nm*Sx2-(Sx)^2)-Sx3*(nm*Sx3-Sx*Sx2)+Sx2*(Sx3*Sx-(Sx)^2));

ea2 = (Sx2y*(nm*Sx2-(Sx)^2)-Sx3*(nm*Sxy-Sx*Sy)+Sx2*(Sxy*Sx-Sy*(Sx)^2))/aux;
ea1 = (Sx4*(nm*Sxy-Sx*Sy)-Sx2y*(nm*Sx3-Sx2*Sx)+Sx2*(Sx3*Sy-Sx2*Sxy))/aux;
ea0 = (Sx4*(Sx2*Sy-Sx*Sxy)-Sx3*(Sx3*Sy-Sx2*Sxy)+Sx2y*(Sx3*Sx-(Sx2)^2))/aux;
ey = ea2*x.^2+ea1*x+ea0;

epy1 = a2*px.^2+a1*px+a0;
epy2 = ea2*px.^2+ea1*px+ea0;
er1 = epy1-py;
MSE1 = sqrt(er1'*er1)/nm;
er2 = epy2-py;
MSE2 = sqrt(er2'*er2)/nm;

figure(1);
plot(x,y, 'r-');
grid on;
hold on;
plot(px, py, 'bo');
plot(x, ey, 'g-');