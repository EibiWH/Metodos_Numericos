clc;
close all;
clear all;

xi = -10;
xf = 10;
dx = 0.001;
er = 0.2;

a1 = 1.57;
a0 = -3;

x = xi:dx:xf;
y = a1*x+a0;
rg = y(end)-y(1);

nm = 20;
px = sort((xf - xi)*rand(nm,1)+xi);
py = a1*px+a0+2*rg*er*rand(nm,1)-rg*er;

Sx = ones(1,nm)*px;
Sy = ones(1,nm)*py;
Sxy = px'*py;
Sx2 = px'*px;

ea1 = (nm*Sxy-Sx*Sy)/(nm*Sx2-Sx2);
ea0 = Sy/nm-ea1*Sx/nm;
ey = ea1*x+ea0;

epy1 = a1*px+a0;
epy2 = ea1*px+a0;
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