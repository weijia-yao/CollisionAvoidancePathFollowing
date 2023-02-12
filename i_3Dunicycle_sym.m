%% symbolic expression
clc;
E = [0, -1; 1, 0];

syms x y z;
k11 = 2; k12 = 2;
e11 = 0.035*log(((x + 3.0)^2 + y^2)^(1/2) + 1.0)*((x + 3.0)^2 + y^2) - 0.048*log(((x - 1.5)^2 + y^2)^(1/2) + 1.0)*((x - 1.5)^2 + y^2) - 0.048*log(((y - 1.3)^2 + (x + 0.75)^2)^(1/2) + 1.0)*((y - 1.3)^2 + (x + 0.75)^2) - 0.048*log(((y + 1.3)^2 + (x + 0.75)^2)^(1/2) + 1.0)*((y + 1.3)^2 + (x + 0.75)^2) + 0.035*log(((y - 2.6)^2 + (x - 1.5)^2)^(1/2) + 1.0)*((y - 2.6)^2 + (x - 1.5)^2) + 0.035*log(((y + 2.6)^2 + (x - 1.5)^2)^(1/2) + 1.0)*((y + 2.6)^2 + (x - 1.5)^2) - 1.0
e12 = z
n11 = [diff(e11,x); diff(e11,y); diff(e11,z)]
n12 = [diff(e12,x); diff(e12,y); diff(e12,z)]
v1 = cross(n11, n12) - k11*e11*n11 - k12*e12*n12;

k2 = 1;
e2 = (x+2.8)^2 + y^2 +z^2 - 1
n2 = [diff(e2,x); diff(e2,y); diff(e2,z)]
v2 = cross(n2,[1;0;0]) - k2*e2*n2;


c = -0.72;
l1 = 0.1; l2 = 0.1;
f1 = exp( l1/(c - e2) );
f2 = exp( l2/(e2 - 0) );
bump1 = f1/(f1+f2);
bump2 = 1-bump1;
vfcom = bump1*v1/norm(v1) + bump2*v2/norm(v2);

diary on
Jv2 = [diff(v2, x) diff(v2, y) diff(v2,z)]
Jvfcom = [diff(vfcom, x) diff(vfcom, y) diff(vfcom,z)]
Jv1 = [diff(v1, x) diff(v1, y) diff(v1,z)]
diary off




