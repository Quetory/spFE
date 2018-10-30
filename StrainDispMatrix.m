function [B,dJ]=StrainDispMatrix(DN,XYZ)
NPE = size(XYZ,1);

J = DN*XYZ;
dJ = det(J);

DNj = J\DN;

if NPE == 4
    
    B(1,1:2:2*NPE) = DNj(1,:);
    B(2,2:2:2*NPE) = DNj(2,:);
    
    B(3,1:2:2*NPE) = DNj(2,:);
    B(3,2:2:2*NPE) = DNj(1,:);
    
elseif NPE == 8
    B(1,1:3:3*NPE) = DNj(1,:);
    B(2,2:3:3*NPE) = DNj(2,:);
    B(3,3:3:3*NPE) = DNj(3,:);

    B(4,1:3:3*NPE) = DNj(2,:);
    B(4,2:3:3*NPE) = DNj(1,:);

    B(5,2:3:3*NPE) = DNj(3,:);
    B(5,3:3:3*NPE) = DNj(2,:);

    B(6,1:3:3*NPE) = DNj(3,:);
    B(6,3:3:3*NPE) = DNj(1,:);

elseif NPE==20      
    
    B(1,1:3:3*NPE) = DNj(1,:);
    B(2,2:3:3*NPE) = DNj(2,:);
    B(3,3:3:3*NPE) = DNj(3,:);

    B(4,1:3:3*NPE) = DNj(2,:);
    B(4,2:3:3*NPE) = DNj(1,:);

    B(5,2:3:3*NPE) = DNj(3,:);
    B(5,3:3:3*NPE) = DNj(2,:);

    B(6,1:3:3*NPE) = DNj(3,:);
    B(6,3:3:3*NPE) = DNj(1,:);

end

