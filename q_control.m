function [K,K0]=q_control(G,Q)

G=ss(G);
[n,m]=size(G.b);
[p,n]=size(G.c);
A=G.a;
B=G.b;
C=G.c;
D=G.d;

Qweight=G.c'*G.c;
Rweight=eye(m);
F=-lqr(A,B,Qweight,Rweight);
L=-lqr(A',C',Qweight,Rweight)';

K0=ss(A+B*F+L*C+L*D*F,[-L B+L*D],[F;-(C+D*F)],...
    [zeros(m,p) eye(m);eye(p) -D]);

K=lft(K0,Q,p,m);