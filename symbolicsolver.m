syms t1 t2 r1 r2 qdot k

A = [ln(r1) 1;
     ln(r2) 1];
B = [t1 - (qdot*r1^2)/(4*k); T2 - (qdot*r2^2)/(4*k)];

X = linsolve(A,B)