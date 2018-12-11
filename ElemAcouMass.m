function [Melem] = ElemAcouMass(XYZ,mat)

if ~isempty(mat.cs)
    Keff = mat.K_s;
    rho0= mat.rho_s;
else
    rho0 = mat.rho;
    Keff  = mat.c^2*mat.rho;
end

NPE = size(XYZ,1);

if NPE == 8
    IP.XI = [-1  -1  -1;...
              1  -1  -1;...
              1   1  -1;...
             -1   1  -1 ;...
             -1  -1   1;...
              1  -1   1;...
              1   1   1;...
             -1   1   1]/sqrt(3);
    IP.WT = ones(1,8);
elseif NPE == 20
    IP.WT  = [ 0.335180055401662*ones(1,8) 0.886426592797784*ones(1,6)];
    p1 = 0.758786910639328;
    p2 = 0.795822425754222;
    p3 = 0;
    IP.XI = [  -p1  p1  p1 -p1 -p1  p1  p1 -p1  p3  p3 p2  p3 -p2  p3;... 
                -p1 -p1  p1  p1 -p1 -p1  p1  p1  p3 -p2 p3  p2  p3  p3;...
                -p1 -p1 -p1 -p1  p1  p1  p1  p1 -p2  p3 p3  p3  p3  p2].';

end

Melem = zeros(NPE,NPE);


for IN = 1:numel(IP.WT)
    XI = IP.XI(IN,:);
    WT = IP.WT(IN);
    if NPE == 8
        [N,DN] = shape3D8(XI);
    elseif NPE == 20
        [N,DN] = shape3D20(XI);
    end
    J    = DN*XYZ;
    dJ   = det(J);
    
    Melem = Melem + 1/Keff*WT*(N.')*(N)*dJ;
end
Melem = rho0*Melem;