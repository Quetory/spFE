clear,close all
clc

XYZ = [-1  -1  -1;...
        1  -1  -1;...
        1   1  -1;...
       -1   1  -1 ;...
       -1  -1   1;...
        1  -1   1;...
        1   1   1;...
       -1   1   1;...
        0  -1  -1;...
        1   0  -1;...
        0   1  -1;...
       -1   0  -1;...
        0  -1   1;...
        1   0   1;...
        0   1   1;...
       -1   0   1;...
       -1  -1   0;...
        1  -1   0;...
        1   1   0;...
       -1   1   0 ];
   
% idx = find(XYZ(:,3)==-1);
          
[Ct] = ElemThermalMass(XYZ,2700,800);

[Kt] = ElemThermalConductivity(XYZ,110);

% [Ktm] = ElemThermalMassflow(XYZ, [1 0 0]);


Ct = Ct + (Ct.'-Ct)/2;
Kt = Kt + (Kt.'-Kt)/2;

% Ct(:,idx)=[];
% Ct(idx,:)=[];
% Kt(:,idx)=[];
% Kt(idx,:)=[];

cik = Ct\Kt;

norm(cik.'-cik)

return
%%

%%
[V,l]=eig(Mt);

l = real(diag(l));

figure(1)
clf
for t = 0:1:100

    scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),[],real(V(:,end))*sin(2*pi/10*t))
    
    xlim([-2 2])
    ylim([-2 2])
    zlim([-2 2])
    caxis([ min(V(:)) max(V(:))])
    drawnow
%     pause
end
    


return
%% Read in Stiffness matrix
% [ K, b0, info.K ] = ansys_read_hb_binary('SM_f_k.hb', 1); 

[ CK, b1, ~ ] = ansys_read_hb_binary('thermalmass_2//SM_f_ck.hb', 1);   % C+K
[ CK2, ~, ~ ] = ansys_read_hb_binary('thermalmass_2//SM_f_2ck.hb', 1); % 2*C+K

% [ ~, b2, ~ ] = ansys_read_hb_binary('SM_f_ck.hb', 1);   % C+K


% Extract System matrices
C = CK2-CK;
K = CK-C;




%%

[w,l]=eig(full(C))

n0 = rank(null(full(C)))
 w.'*C*w
L = diag(l);
Kw = w.'*K*w;
K00 = Kw(1:n0,1:n0);
K0C = Kw(1:n0,n0+1:end);
KC0 = Kw(n0+1,1:n0);
KCC = Kw(n0+1:end,n0+1:end);


Tg = [-K00\K0C; eye(size(C,1)-n0)];

Kwt = Tg.'*Kw*Tg;



