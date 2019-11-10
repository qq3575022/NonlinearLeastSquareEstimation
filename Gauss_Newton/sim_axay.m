clc, clear, close all

% ========================  Load instantaneous trilateration r rdot a  ==========================
load('yVectorData.mat')
h = 0.001/100;

% =================== Load Ground truth coordinates of the moving tag path ===========
load('X.mat')
load('A.mat')

w_a = A(1,1:190991);  W_A = NaN(1,length(tf)); 
w_v = X(2,1:190991);  W_V = NaN(1,length(tf));
w   = X(1,1:190991);  W   = NaN(1,length(tf));

k = 1;
for n = 2:length(w)
   if mod(n,710) == 0
        k      = k + 1;
        W_A(k) = w_a(n);
        W_V(k) = w_v(n);
        W(k)   = w(n);
    end   
end

W_A(1) = w_a(1); W_V(1) = w_v(1); W(1) = w(1);

x     = NaN(3,length(w));     y = NaN(3,length(w));     XX = NaN(6,length(tf)); 
x_num = NaN(3,length(w)); y_num = NaN(3,length(w)); XX_num = NaN(6,length(tf)); 

x(1,1) = -0.3*sin(w(1)) + 1.95; x(2,1) = 0; x(3,1) = 0;
y(1,1) =  0.3*cos(w(1)) + 1.51; y(2,1) = 0; y(3,1) = 0;

ax = NaN(1,length(w));    ay = NaN(1,length(w));

k = 1;
for n = 2:length(w)-1
    x(1,n) =  0.3*sin(-w(n)) + 1.95;%2;
    x(2,n) = -0.3*cos(-w(n))*w_v(n);
    x(3,n) = -0.3*sin(-w(n))*w_v(n)^2 - 0.3*cos(-w(n))*w_a(n);
    
    y(1,n) =  0.3*cos(-w(n)) + 1.51;%1.5;
    y(2,n) =  0.3*sin(-w(n))*w_v(n);
    y(3,n) = -0.3*cos(-w(n))*w_v(n)^2 + 0.3*sin(-w(n))*w_a(n);
    
    ax(1,n) =  -x(3,n)*sin(w(n))+y(3,n)*cos (w(n));
    ay(1,n) =  -x(3,n)*cos(w(n))-y(3,n)*sin(w(n));
    
    x_num(1,n) = x(1,n);
    x_num(2,n) = (x(1,n) - x(1,n-1))/h;
    x_num(3,n) = (x(2,n) - x(2,n-1))/h;
    
    y_num(1,n) = y(1,n);
    y_num(2,n) = (y(1,n) - y(1,n-1))/h;
    y_num(3,n) = (y(2,n) - y(2,n-1))/h;
    
    
    if mod(n,710) == 0
        k = k+1;  
        
        XX(1,k) = x(1,n); 
        XX(2,k) = x(2,n); 
        XX(3,k) = x(3,n);
        
        XX(4,k) = y(1,n); 
        XX(5,k) = y(2,n); 
        XX(6,k) = y(3,n);
        
        XX_num(1,k) = x_num(1,n); 
        XX_num(2,k) = x_num(2,n); 
        XX_num(3,k) = x_num(3,n);
        
        XX_num(4,k) = y_num(1,n); 
        XX_num(5,k) = y_num(2,n); 
        XX_num(6,k) = y_num(3,n);
        
        AXY(1,k) = ax(1,n);
        AXY(2,k) = ay(1,n);
    end
end

XX(1,1) = x(1,1);     XX(2,1) = 0;     XX(3,1) = 0;
XX(4,1) = y(1,1);     XX(5,1) = 0;     XX(6,1) = 0; 

XX_num(1,1) = x(1,1); XX_num(2,1) = 0; XX_num(3,1) = 0;
XX_num(4,1) = y(1,1); XX_num(5,1) = 0; XX_num(6,1) = 0; 

figure
subplot(3,1,1);plot(tf,W_A);title('Grount Truth Angular Acceleration')
subplot(3,1,2);plot(tf,W_V);title('Grount Truth Angular Velocity')
subplot(3,1,3);plot(tf,W);  title('Grount Truth Angular Position')

figure
subplot(6,1,1);plot(tf,XX(1,:));title('X Position Coordinates Calculated from Formula')
subplot(6,1,2);plot(tf,XX(2,:));title('X Velocity Calculated from Formula')
subplot(6,1,3);plot(tf,XX(3,:));title('X Acceleration Calculated from Formula')
subplot(6,1,4);plot(tf,XX(4,:));title('Y Position Coordinates Calculated from Formula')
subplot(6,1,5);plot(tf,XX(5,:));title('Y Velocity Calculated from Formula')
subplot(6,1,6);plot(tf,XX(6,:));title('Y Acceleration Calculated from Formula')

figure
subplot(6,1,1);plot(tf,XX_num(1,:));
subplot(6,1,2);plot(tf,XX_num(2,:));title('X Velocity Calculated Numerically')
subplot(6,1,3);plot(tf,XX_num(3,:));title('X Acceleration Calculated Numerically')

subplot(6,1,4);plot(tf,XX_num(1,:));
subplot(6,1,5);plot(tf,XX_num(5,:));title('Y Velocity Calculated Numerically')
subplot(6,1,6);plot(tf,XX_num(6,:));title('Y Acceleration Calculated Numerically')

figure
subplot(2,1,1);plot(tf,AXY(1,:),tf,ax1);legend('Calculated Based on Ground Truth Angle', 'Measurement');title('x^B Acceleration')
subplot(2,1,2);plot(tf,AXY(2,:),tf,ay1);legend('Calculated Based on Ground Truth Angle', 'Measurement');title('y^B Acceleration')