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

nm = 15;
no = 1;
px = sort((xf-xi)*rand(nm,1)+xi);
py = a1*px+a0+2*rg*er*rand(nm,1)-rg*er;
id = randi(nm, no, 1);
po = 4*rg*rand(no,1)-2*rg;
po(po<py(end)&po>py(1)) = po(po<py(end)&po>py(1))+(2*(po(po<py(end)&po>py(1))>0)-1)*rg;
py(id) = po;

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

err2 = er2.*er2;
[merr2, ide]=max(err2);
 
pxo = [px(1:ide-1); px(ide+1:nm)];
pyo = [py(1:ide-1); py(ide+1:nm)];
nmo = nm-1;
 
Sxo = ones(1,nmo)*pxo;
Syo = ones(1,nmo)*pyo;
Sxyo = pxo'*pyo;
Sx2o = pxo'*pxo;
 
ea1o = (nmo*Sxyo-Sxo*Syo)/(nmo*Sx2o-Sxo^2);
ea0o = Syo/nmo-ea1o*Sxo/nmo;
eyo = ea1o*x+ea0o;
 
epy3 = ea1o*pxo+ea0o;
er3 = epy3-pyo;
MSE3 = sqrt(er3'*er3)/nmo;

figure(1);
plot(x, y, 'r-');
grid on;
hold on;
plot(px, py, 'bo');
plot(x, ey, 'g-');
plot(px(ide), py(ide), 'k*')
plot(x, eyo, 'k-');