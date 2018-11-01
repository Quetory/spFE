function [Kelem] = ElemStiffness(XYZ, E, nu,eas)
%% Input Parsing
if isempty(eas)
    eas = 0;
end

%% Define Gauss integration points
NPE = size(XYZ,1);
NDOF = size(XYZ,2);

if NPE==4
    IP.XI = [-1  -1 ;...
              1  -1 ;...
              1   1 ;...
             -1   1  ]/sqrt(3);
    IP.WT = ones(1,4);

elseif NPE == 8 
    IP.XI = [-1  -1  -1;...
              1  -1  -1;...
              1   1  -1;...
             -1   1  -1 ;...
             -1  -1   1;...
              1  -1   1;...
              1   1   1;...
             -1   1   1]/sqrt(3);
    IP.WT = ones(1,8);

elseif NPE == 20 || eas
    IP.WT  = [ 0.335180055401662*ones(1,8) 0.886426592797784*ones(1,6)];
    p1 = 0.758786910639328;
    p2 = 0.795822425754222;
    p3 = 0;
    IP.XI = [  -p1  p1  p1 -p1 -p1  p1  p1 -p1  p3  p3 p2  p3 -p2  p3;... 
               -p1 -p1  p1  p1 -p1 -p1  p1  p1  p3 -p2 p3  p2  p3  p3;...
               -p1 -p1 -p1 -p1  p1  p1  p1  p1 -p2  p3 p3  p3  p3  p2].';

end


if NDOF==2 
     D = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];                 %assuming plane stress
     D = E/((1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 1-2*nu];    %assuming plane strain
elseif NDOF==3
    D = diag( [1-nu, 1-nu,1-nu, 1/2-nu, 1/2-nu, 1/2-nu])/2;
    D = D + diag([nu,nu,0,0,0],1) + diag([nu,0,0,0],2);
    D= D.'+D;
    D = E/((1+nu)*(1-2*nu))*D;
end
   
        
if eas
    [T0, dJ0] = EAS_M_transformation(XYZ,[0 0 0]);
    T0 = T0*dJ0;
    m(2) = 4;
    m(3) = 11;
%     m(3) = 17;
    L = zeros(m(NDOF),NDOF*NPE);
    Q = zeros(m(NDOF));
else
    T = [];
    L = [];
    Q = [];
    
end
%% Numerical integration to get Kelem
Kelem = zeros(NDOF*NPE,NDOF*NPE);

for IN = 1:numel(IP.WT)
    XI = IP.XI(IN,:);
    WT = IP.WT(IN);
    if NPE==4 && ~eas
        [~,DN] = shape2D4(XI);
        [B,dJ] = StrainDispMatrix(DN,XYZ); 
    elseif NPE ==4 && eas
        [~,DN] = shape2D4(XI);
        [B,dJ] = StrainDispMatrix(DN,XYZ); 
        [M]    = EAS_Matrix(XI,T0,dJ);
        L = L + WT*(M.')*D*B*dJ;
        Q = Q + WT*(M.')*D*M*dJ;
        
    elseif NPE == 8 && ~eas
        [~,DN] = shape3D8(XI);
        [B,dJ] = StrainDispMatrix(DN,XYZ);

        
    elseif NPE == 8 && eas
        [~,DN] = shape3D8(XI);
%         [B,dJ] = StrainDispMatrix(DN,XYZ);
        [Bans,dJ] = ANS_interpolation(DN,XI,XYZ);
        [Ti,~] = EAS_M_transformation(XYZ,XI);
        B= Ti*Bans;
        
        [M]    = EAS_Matrix(XI,T0,dJ);
        
        L = L + WT*(M.')*D*B*dJ;
        Q = Q + WT*(M.')*D*M*dJ;
        
    elseif NPE == 20
        [~,DN] = shape3D20(XI);
        [B,dJ] = StrainDispMatrix(DN,XYZ);
        
    end
   
   Kelem = Kelem + WT*B.'*D*B*dJ ;
   
end

if eas % Static condensation step for EAS 
    Kelem = Kelem - L.'*(Q\L);
end   