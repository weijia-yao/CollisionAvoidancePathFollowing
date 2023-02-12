% symbolic expression
syms x y;
R = 2; a = 1; b = 0.5; E = [0, -1; 1, 0];
k1 = 1; k2 = 1;

phi = x^2 + y^2 - R^2;
varphi = ((x-0.9)^2+(y-R)^2)*((x+0.9)^2+(y-R)^2) - 0.9;

nphi    = [diff(phi,x); diff(phi,y)];
nvarphi = [diff(varphi,x); diff(varphi,y)]

v1 = E * nphi - k1 * phi * nphi;
v2 = E * nvarphi - k2 * varphi * nvarphi;

%% bump functions
c = -0.2;
kk = 0.1;
f1 = piecewise(varphi <= c, 0, varphi>c, exp(kk/(c-varphi)));
f2 = piecewise(varphi < 0, exp(kk*(1 + 1/varphi)), varphi >= 0, 0);
bump1 = f1/(f1+f2);
bump2 = f2/(f1+f2);

figure; xlabel('x'); ylabel('y'); zlabel('z'); % axis equal;
fs1 = fsurf(bump1, [-1.1 1.1 1.3 2.7]);
hold on; 
fs2 = fsurf(bump2, [-1.1 1.1 1.3 2.7]);  
fs2.FaceAlpha=0.4; 
title('bump function 1 and 2')
hold off;
figure; axis equal; xlabel('x'); ylabel('y'); zlabel('z');
fs3 = fsurf(bump1+bump2, [-1.1 1.1 1.3 2.7]); 
title('Sum of bump function 1 and 2')