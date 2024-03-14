%
% [SM,U,V]=smform(G)
% 
% finding the smith mcmillan form of a MIMO LTI Object G
%
% M=U*G*V
% 
function [M,U,V,Gsym,Nsym,Dsym]=smform(G)

% put it in state space form
G=ss(G);
% convert to symbolic form
s=sym('s');
n=order(G);
% symbolic numerator and denominator
Nsym=G.C*adjoint(s*eye(n)-G.A)*G.B;
Dsym=det(s*eye(n)-G.A);
Gsym=simplify(Nsym/Dsym);
% find Smith form of N
[U,V,Nsym_smith]=smithForm(Nsym,s);
% divide by d(s) and simplify
M=simplify(Nsym_smith/Dsym);


