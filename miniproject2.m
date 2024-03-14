close all

c1=.01; c2=.005; c3=.001; l1=.001; l2=.05; l3=.001; r1=100; r2=100;
A = [0 0 0 1/c1 0 0;
     0 -1/r2/c2 0 1/c2 -1/c2 -1/c2;
     0 0 0 0 0 1/c3;
     -1/l1 -1/l1 0 0 0 0;
     0 1/l2 0 0 -1/r2/l2 0;
     0 1/l3 -1/l3 0 0 0];
B2 = [0 0; 1/r1/c2 0; 0 0; 1/l1 0;0 1/l2; 0 0];
C2 = [0 0 0 1 0 0; 0 1/r1 0 0 0 0];
D22=[0 0; -1/r1 0];
G_ss = ss(A,B2,C2,D22);

%Full P matrix
e1 = 1; e2 = 1;
B1 = [0 4 0 0;
    0 0 0 0; 
    0 0 0 0;
    0 0 0 0;
    0 0 0 0;
    0 0 0 0];
C1 = [G_ss.C(1,:);0 0 0 0 0 0;0 0 0 0 0 0];
D11 = [-1 0 0 0;
       0 0 0 0;
       0 0 0 0];
D12 = [0 0; 
       e1 0;
       0 e2];
D21 = [-1 0 1 0;
       0 0 0 1];
B = cat(2,B1,B2);
C = cat(1,C1,C2);
D = cat(1,cat(2,D11,D12),cat(2,D21,D22));

P = ss(A,B,C,D)

% Smith Mcmillan

[M,U,V,Gsym,Nsym,Dsym] = smform(G_ss);

% Coprime factorization
epsilon = sym(zeros(2,2));
psi = sym(zeros(2,2));

for i = 1:2
    [n,d] = numden(M(i,i));
    epsilon(i,i) = n;
    psi(i,i) = d;
end

Nr = simplify(inv(U)*epsilon);
Dr = simplify(V*psi);
Nl = simplify(epsilon*inv(V));
Dl = simplify(psi*U);

right_check = simplify(Nr*inv(Dr)-Gsym);
left_check = simplify(inv(Dl)*Nl-Gsym);

%stable coprime factorization

[factr,Mr,Nr] = rncf(G_ss);
[factl,Ml,Nl] = lncf(G_ss);

% Controllability and Observiabilty gramiens to find best input/output pair

controllability = [eig(gram(G_ss(1,1),'c')),eig(gram(G_ss(1,2),'c')),eig(gram(G_ss(2,1),'c')),eig(gram(G_ss(2,2),'c'))];

observability = [eig(gram(G_ss(1,1),'o')),eig(gram(G_ss(1,2),'o')),eig(gram(G_ss(2,1),'o')),eig(gram(G_ss(2,2),'o'))];

poles = pole(G_ss);
tzeros = tzero(G_ss);


%Parameterization

E = [0 0 0 c1 0 0;
     0 -c2/r1 0 c2 -c2 -c2;
     0 0 0 0 0 c3;
     -l1 -l1 0 0 0 0;
     0 l2 0 0 -l2/r2 0;
     0 l3 -l3 0 0 0];
Bdelbar = [0 0;0 1;0 0;1 0;0 0;0 0];
Cdel = [-1 -1 0 0 0 0;0 -1 0 1 -1 -1];
Bdel = inv(E)*Bdelbar;
Ddeldel = [1 0;0 1];
Ddel1 = [-1 0 1 0;0 -1/r1 0 1];
Ddel2 = [0 0;-1/r1 0];
D1del = zeros(3,2);
D2del= zeros(2);
B = cat(2,Bdel,B);
C = cat(1,Cdel,C);
D = cat(2,cat(1,cat(1,Ddeldel,D1del),D2del),cat(1,cat(2,Ddel1,Ddel2),D));

P = ss(A,B,C,D);


L1 = ureal('l1',.001,'Percentage',30);
C2 = ureal('c2',.005,'Percentage',30);

delta = [L1,0;0,C2];

P_uncertain = lft(delta,P)
Pnomval = P_uncertain.NominalValue

%Stabalizing Q control
[KQ,K0] = q_control(G_ss,ones(2,2));
figure,bode(KQ)

%nyquist plots
w=logspace(-3,3,100);
F = feedback(KQ*G_ss,eye(2));
gamma = mvar_nyquist(KQ*G_ss,w);
figure(1); 
plot(real(gamma(1,:)),imag(gamma(1,:)),'b',...
    real(gamma(1,:)),-imag(gamma(1,:)),'b:',...
    real(gamma(2,:)),imag(gamma(2,:)),'r',...
    real(gamma(2,:)),-imag(gamma(2,:)),'r:',...
    'linewidth',2);
grid
figure,nyquist(KQ*G_ss)

%Passivity therorem
passivity = norm((eye(2)-G_ss) * inv(eye(2)+G_ss),'inf')

