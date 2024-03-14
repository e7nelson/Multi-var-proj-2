%
% calculating the left and right stable coprime factorization 
% of a proper rational matrix G(s)
% 

function [M,N,Y,X,Mt,Nt,Xt,Yt]=scf(G)

G=ss(G);
[n,m]=size(G.b);
[p,n]=size(G.c);
A=G.a;
B=G.b;
C=G.c;
D=G.d;

Q=G.c'*G.c;
R=eye(m);
F=-lqr(A,B,Q,R);
L=-lqr(A',C',Q,R)';
M=ss(A+B*F,B,F,eye(m));
N=ss(A+B*F,B,C+D*F,D);
Mt=ss(A+L*C,L,C,eye(p));
Nt=ss(A+L*C,B+L*D,C,D);
X=ss(A+B*F,-L,C+D*F,eye(p));
Y=ss(A+B*F,-L,F,zeros(m,p));
Xt=ss(A+L*C,-(B+L*D),F,eye(m));
Yt=ss(A+L*C,L,F,zeros(m,p));
end