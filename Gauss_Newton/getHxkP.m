function H = getHxkP(x)
x1 = [-0.05, 1.5];
x2 = [2, 3.0];
x3 = [2.7, 0.05];
H = NaN(3,1);

%x: x(1); xdot: x(2); y: x(4); ydot: x(5)
H(1,1) = sqrt((x(1)-x1(1))^2+(x(4)-x1(2))^2);

H(2,1) = sqrt((x(1)-x2(1))^2+(x(4)-x2(2))^2);

H(3,1) = sqrt((x(1)-x3(1))^2+(x(4)-x3(2))^2);
end