function H = getHxkPVA(x)
x0 = [1.95, 1.51];

x1 = [-0.05, 1.5]; 
x2 = [2, 3.0]; 
x3 = [2.7, 0.05];

H = NaN(10,1);

H(1,1) = sqrt((x(1)-x1(1))^2+(x(4)-x1(2))^2);
H(2,1) = ((x(1)-x1(1))*x(2)+(x(4)-x1(2))*x(5))/sqrt((x(1)-x1(1))^2+(x(4)-x1(2))^2);

H(3,1) = sqrt((x(1)-x2(1))^2+(x(4)-x2(2))^2);
H(4,1) = ((x(1)-x2(1))*x(2)+(x(4)-x2(2))*x(5))/sqrt((x(1)-x2(1))^2+(x(4)-x2(2))^2);

H(5,1) = sqrt((x(1)-x3(1))^2+(x(4)-x3(2))^2);
H(6,1) = ((x(1)-x3(1))*x(2)+(x(4)-x3(2))*x(5))/sqrt((x(1)-x3(1))^2+(x(4)-x3(2))^2);

H(7,1) =  x(7);%wrapToPi(x(7));
H(8,1) =  x(8);

H(9,1) = -x(3)*sin(x(7)) + x(6)*cos(x(7));
H(10,1)= -x(3)*cos(x(7)) - x(6)*sin(x(7));

% H(7,1) =  x(3);%x(3)*x(4)/sqrt(x(1)^2+x(4)^2)+x(6)*x(1)/sqrt(x(1)^2+x(4)^2);
% H(8,1) =  x(6);%-x(3)*x(1)/sqrt(x(1)^2+x(4)^2)+x(6)*x(4)/sqrt(x(1)^2+x(4)^2);

end