clear,close all
clc

%%
XY = [-1 -1 ;...
        1 -1 ;...
        1 1 ;...
       -1 1 ];
     
E=1;
nu=0.3;
eas=0;
rho = 1;

[K,~] = ElemStiffness(XY, E, nu,eas);


[Keas,~] = ElemStiffness(XY, E, nu,1);
[Melem] = ElemMass(XY,rho);

[V,L]=eig(Keas,Melem);
[V1,L1]=eig(K,Melem);

%% 
i= 4
dx = reshape(1/5*V(:,i).' , 2,4).';%

%%
figure

plot(XY(:,1),XY(:,2),'ko' )
grid on
xlim([-2,2])
ylim([-2,2]);

hold all
plot(XY(:,1)+dx(:,1),XY(:,2)+dx(:,2),'r*' )
title(num2str(L(i,i)))