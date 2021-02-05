function [F0,dJ] = EAS_M_transformation(XYZ,XI, ab)
NPE = size(XYZ,1);

%reference: The Finite Element Method: Its Basis and Fundamentals 6th, p373
 
if NPE==8
    if ab==1
        a = 1;
        b=2;
    else
        a=2;
        b=1;
    end
    [~,DN] = shape3D8(XI);
    J = DN*XYZ;
    Ji = (eye(3)/J) ;
    dJ = det(J);

%     T = J0i.';
    T = Ji;

    %     e = (e11 e22 e33 e12 e23 e13 )          
    T11 = (T.*T).';
    
    T12 = a*[ T(1,1)*T(2,1) T(2,1)*T(3,1) T(1,1)*T(3,1);...
              T(1,2)*T(2,2) T(2,2)*T(3,2) T(1,2)*T(3,2);...
              T(1,3)*T(2,3) T(2,3)*T(3,3) T(1,3)*T(3,3)];
            
    T21 =  b*[ T(1,1)*T(1,2) T(2,1)*T(2,2) T(3,1)*T(3,2);...
              T(1,2)*T(1,3) T(2,2)*T(2,3) T(3,2)*T(3,3);...
              T(1,1)*T(1,3) T(2,1)*T(2,3) T(3,1)*T(3,3)];

    T22 =   [ T(1,1)*T(2,2)+T(1,2)*T(2,1) T(2,1)*T(3,2)+T(2,2)*T(3,1) T(1,1)*T(3,2)+T(1,2)*T(3,1)  ;...
              T(1,2)*T(2,3)+T(1,3)*T(2,2) T(2,2)*T(3,3)+T(2,3)*T(3,2) T(1,2)*T(3,3)+T(1,3)*T(3,2)   ;...
              T(1,1)*T(2,3)+T(1,3)*T(2,1) T(2,1)*T(3,3)+T(2,3)*T(3,1) T(1,1)*T(3,3)+T(1,3)*T(3,1)  ];

%     e = (e11 e22 e33 e12 e13 e23 )          
%     T12 = 2*[ T(1,1)*T(2,1) T(1,1)*T(3,1) T(2,1)*T(3,1);...
%               T(1,2)*T(2,2) T(1,2)*T(3,2) T(2,2)*T(3,2);...
%               T(1,3)*T(2,3) T(1,3)*T(3,3) T(2,3)*T(3,3)];
%                        
%     T21 =   [ T(1,1)*T(1,2) T(2,1)*T(2,2) T(3,1)*T(3,2);...
%               T(1,1)*T(1,3) T(2,1)*T(2,3) T(3,1)*T(3,3);...
%               T(1,2)*T(1,3) T(2,2)*T(2,3) T(3,2)*T(3,3)];
%          
    F0 = [T11 T12; T21 T22];
elseif NPE==4
    %reference: The Finite Element Method: Its Basis and Fundamentals 6th, p373
     [~,DN] = shape2D4(XI);
     J0 = DN*XYZ;
     J0i = (eye(2)/J0) ;
     dJ0 = det(J0);
     T = J0i.';
     
     F0 = [             (T.*T).'          [2*T(1,1)*T(2,1); 2*T(1,2)*T(2,2)]  ;...
                T(1,1)*T(1,2) T(2,1)*T(2,2)    T(1,1)*T(2,2)+T(1,2)*T(2,1)     ];
     
end
        
     
     
     