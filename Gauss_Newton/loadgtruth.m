function [W_A, W_V, W, XX, rr, AXY, Orient] = loadgtruth(td)
%The x and y coordinates of tag
x1 = [-0.05, 1.5]; x2 = [2, 3.0]; x3 = [2.7, 0.05];

load('X.mat');load('A.mat');load('AxB.mat');load('AyB.mat')

w_a = A(1,1:190991);     w_v = X(2,1:190991);     w = X(1,1:190991);     orit = X(1,1:190991);  
W_A = NaN(1,length(td)); W_V = NaN(1,length(td)); W = NaN(1,length(td)); Orient = NaN(1,length(td));

k = 1;
for n = 2:length(w)
   if mod(n,710) == 0
        k      = k + 1;
        W_A(k) = w_a(n);
        W_V(k) = w_v(n);
        W(k)   = w(n);
        
        Orient(k) = orit(n);
    end   
end

offset = 0;
for n = 2:1:length(Orient)
    if Orient(n) - Orient(n-1) < -2
        offset = 2*pi;
    end
    Orient(n) = Orient(n) + offset;
        
end

W_A(1) = w_a(1); W_V(1) = w_v(1); W(1) = w(1);

x = NaN(3,length(w)); y = NaN(3,length(w)); XX = NaN(6,length(td)); 
x(1,1) = -0.3*sin(w(1)) + 1.95; x(2,1) = 0; x(3,1) = 0;
y(1,1) =  0.3*cos(w(1)) + 1.51; y(2,1) = 0; y(3,1) = 0;

ax = NaN(1,length(w));ay = NaN(1,length(w));

k = 1;
for n = 2:length(w)-1
    x(1,n) = -0.3*sin(w(n)) + 1.95;%2;
    x(2,n) = -0.3*cos(w(n))*w_v(n);
    x(3,n) =  0.3*sin(w(n))*w_v(n)^2 - 0.3*cos(w(n))*w_a(n);
    
    y(1,n) =  0.3*cos(w(n)) + 1.51;%1.5;
    y(2,n) = -0.3*sin(w(n))*w_v(n);
    y(3,n) = -0.3*cos(w(n))*w_v(n)^2 - 0.3*sin(w(n))*w_a(n);
    
    ax(1,n) = -x(3,n)*sin(w(n))+y(3,n)*cos(w(n));
    ay(1,n) = -x(3,n)*cos(w(n))-y(3,n)*sin(w(n));
%     
    if mod(n,710) == 0
        k = k+1;      
        XX(1,k) = x(1,n); 
        XX(2,k) = x(2,n); 
        XX(3,k) = x(3,n);
        
        XX(4,k) = y(1,n); 
        XX(5,k) = y(2,n); 
        XX(6,k) = y(3,n);
        
        AXY(1,k) = ax(1,n);%AxB(1,n);
        AXY(2,k) = ay(1,n);%AyB(1,n);
    end
end

XX(1,1) = x(1,1); XX(2,1) = 0; XX(3,1) = 0;
XX(4,1) = y(1,1); XX(5,1) = 0; XX(6,1) = 0; 

% get r1 r2 r3
r1 = NaN(2,length(w));  r2  = NaN(2,length(w)); r3  = NaN(2,length(w));
rr = NaN(8,length(td)); axx = NaN(1,length(w)); ayy = NaN(1,length(w));

r1(1,1) = sqrt((x(1,1)-x1(1))^2 + (y(1,1)-x1(2))^2);  r1(2,1) = 0;
r2(1,1) = sqrt((x(1,1)-x2(1))^2 + (y(1,1)-x2(2))^2);  r2(2,1) = 0;
r3(1,1) = sqrt((x(1,1)-x3(1))^2 + (y(1,1)-x3(2))^2);  r3(2,1) = 0;

k = 1;
for n = 2:length(w)
    
    H = getHxkPVA([x(1,n),x(2,n),x(3,n),y(1,n),y(2,n),y(3,n),w(n),w_v(n)]);
    r1(1,n) = H(1);
    r1(2,n) = H(2);
    r2(1,n) = H(3);
    r2(2,n) = H(4);
    r3(1,n) = H(5);
    r3(2,n) = H(6);
    wx(n)   = H(7); % orientation
    wy(n)   = H(8); % angular velocity
    axx(n)  = H(9); % ax
    ayy(n)  = H(10);% ay
    
    if mod(n,710) == 0
        k = k+1;
        rr(1,k) = r1(1,n);rr(2,k) = r1(2,n);
        rr(3,k) = r2(1,n);rr(4,k) = r2(2,n);
        rr(5,k) = r3(1,n);rr(6,k) = r3(2,n);
        rr(7,k) = axx(n); rr(8,k) = ayy(n);
    end
end

rr(1,1) = r1(1,1);rr(2,1) = r1(2,1);
rr(3,1) = r2(1,1);rr(4,1) = r2(2,1);
rr(5,1) = r3(1,1);rr(6,1) = r3(2,1);
rr(7,1) = 0;rr(8,1) = 0;
end