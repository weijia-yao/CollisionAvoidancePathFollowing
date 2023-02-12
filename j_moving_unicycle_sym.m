%% symbolic expression
clc;
E = [0, -1; 1, 0];

syms k1 k2 x y c t obs_v t a b l1 l2;
phi = y - sin(x);

varphi = (x+5-obs_v*t)^2/a^2 + y^2/b^2 - 1;

e1 = phi;
e2 = varphi;
n1 = [diff(phi, x); diff(phi, y)]; 
n2 = [diff(varphi, x); diff(varphi, y)]; 

vf1 = E*n1-k1*e1*n1;
vf2 = E*n2-k2*e2*n2;

f1 = exp( l1/(c - e2) );
f2 = exp( l2/(e2 - 0) );
bump1 = f1/(f1+f2);
bump2 = 1-bump1;
vfcom = bump1*vf1/norm(vf1) + bump2*vf2/norm(vf2);


Jvf2 = [diff(vf2, x) diff(vf2, y)]
Jvfcom = [diff(vfcom, x) diff(vfcom, y)]
Jvf1 = [diff(vf1, x) diff(vf1, y)]


