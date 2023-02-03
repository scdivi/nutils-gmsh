h  = 0.00125;
r1 = 0.1;
r2 = 0.15;
Point(0) = {0,0,0,h};
Point(1) = {r1,0,0,h};
Point(2) = {r2,0,0,h};
Point(3) = {0,r2,0,h};
Point(4) = {0,r1,0,h};
Line(1) = {1,2};
Circle(2) = {2,0,3};
Line(3) = {3,4};
Circle(4) = {4,0,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Physical Line("bottom") = {1};
Physical Line("outer")  = {2};
Physical Line("left")   = {3};
Physical Line("inner")  = {4};
Physical Surface("interior") = {1};
