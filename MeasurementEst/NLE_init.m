clc, clear, close all

% Set timestep
time_stepPos = 0.5;    iternumPos = 40;  
time_step = 1;         iternum    = 20;  

% =============================  Get simulated r rdot a  ============================
global T
h = 1e-4;  tc = 0:h:1.9099;
T = 0.0071;td = 0:T:1.9099;

% Ground truth coordinates of the moving tag path
[W_A, W_V, W, XX, rr, AXY] = getcoordinates(td);

% Parameters of tag
Gt = 14.62;    % tag's antenna gain
X  = 0.85;      % polarization mismatch
M  = 4;         % load modulation factor of the tagoo
f  = 5.8*10^9;

% Parameters of reader
PT = 1;         % reader's transmitted power 1 W
GT = 14.62;     % reader's trasmitter antenna gain 9.5dBi
GR = 14.62;     % reader's receiver   antenna gain 9.5dBi
R  = 15;

% Channel noise error covariance
sigma = 0.0004;%0.0303; 

% Position of readers
x1 = [-0.05, 1.5]; x2 = [2, 3.0]; x3 = [2.7, 0.05];

% Phase cconcatenation
l1 = 0;  l2 = 0;  l3 = 0;

% Get simulated r rdot
z = NaN(2,length(td)-1);       z_prev = NaN(2,length(td)-1); 

% simulated instantaneous trilateration input
z(1,:) = XX(1,2:end);          z_prev(1,:) = XX(1,1:end-1);    % x coordinate
z(2,:) = XX(4,2:end);          z_prev(2,:) = XX(4,1:end-1);    % y coordinate

% simulated instantaneous trilateration r1 r2 r3 r1dot r2dot r3dot
r_sim1 = NaN(1,length(td));    rdot_sim1 = NaN(1,length(td)); 
r_sim2 = NaN(1,length(td));    rdot_sim2 = NaN(1,length(td)); 
r_sim3 = NaN(1,length(td));    rdot_sim3 = NaN(1,length(td)); 

% Concatenate phase and original phase
phi_conc1 = NaN(1,length(td)); phi_mod1 = NaN(1,length(td));
phi_conc2 = NaN(1,length(td)); phi_mod2 = NaN(1,length(td));
phi_conc3 = NaN(1,length(td)); phi_mod3 = NaN(1,length(td));

% Phase difference
diff1 = NaN(1,length(td));diff2 = NaN(1,length(td));diff3 = NaN(1,length(td));

for k = 1:1:length(td)-1
    [phi_conc1(k+1),phi_mod1(k+1),r_sim1(k+1),rdot_sim1(k+1),diff1(k+1),l1] = getsim(x1,f,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,T,l1);
    [phi_conc2(k+1),phi_mod2(k+1),r_sim2(k+1),rdot_sim2(k+1),diff2(k+1),l2] = getsim(x2,f,Gt,M,X,PT,GT,GR,R,sigma,2,k,z,z_prev,T,l2);
    [phi_conc3(k+1),phi_mod3(k+1),r_sim3(k+1),rdot_sim3(k+1),diff3(k+1),l3] = getsim(x3,f,Gt,M,X,PT,GT,GR,R,sigma,3,k,z,z_prev,T,l3); 
end

r_sim1(1) = r_sim1(2);rdot_sim1(1) = 0;phi_mod1(1) = phi_mod1(2);phi_conc1(1) = phi_conc1(2);
r_sim2(1) = r_sim2(2);rdot_sim2(1) = 0;phi_mod2(1) = phi_mod2(2);phi_conc2(1) = phi_conc2(2);
r_sim3(1) = r_sim3(2);rdot_sim3(1) = 0;phi_mod3(1) = phi_mod3(2);phi_conc3(1) = phi_conc3(2);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Verify simulated instantaneous trilateration~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure
subplot(3,1,1),plot(td,r_sim1,'LineWidth',2); hold on; plot(td,rr(1,:),'LineWidth',2);legend('simulated instantaneous trilateration','ground truth');title('simulated r1');
subplot(3,1,2),plot(td,r_sim2,'LineWidth',2); hold on; plot(td,rr(3,:),'LineWidth',2);legend('simulated instantaneous trilateration','ground truth');title('simulated r2');
subplot(3,1,3),plot(td,r_sim3,'LineWidth',2); hold on; plot(td,rr(5,:),'LineWidth',2);legend('simulated instantaneous trilateration','ground truth');title('simulated r3');
xlabel('t [s]')

variancerr = NaN(8,15);
time = NaN(1,25);

[x_meas,A] = rtoxy(r_sim1, r_sim2, r_sim3, AXY,T);

figure
subplot(2,1,1), plot(A(1,:)),title('ax')
subplot(2,1,2), plot(A(2,:)),title('ay')

N = 1;
%==========================  Estimation Based on R:  N = 5===============================

% Simulated simulated instantaneous trilateration y
y = NaN(3*N,length(td)-N);

% Get yN
for l = 1:1:length(td)-N
y(:,l) = getyNP(r_sim1,r_sim2,r_sim3,l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0]*ones(1,length(td)-N);e = NaN(1,length(td)-N);

for m = 1:1:length(td)-N
    
for iter = 1:1:iternumPos
    
E = getEP(y(:,m),x(:,m),N);
orig = norm(E);

JacoN = getJacoN(y(:,m),x(:,m),N,1);
add = pinv(JacoN)*E;

temp = x(:,m) - time_stepPos*add;%x - 0.5.*pinv(JacoN)*E

E_temp = getEP(y(:,m),temp,N);
after = norm(E_temp);

if norm(E_temp) < norm(E)
    x(:,m) = temp;
else
    e(m) = norm(E_temp);
    break
end

end
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0;0,1,T,0,0,0;0,0,1,0,0,0;0,0,0,1,T,0;0,0,0,0,1,T;0,0,0,0,0,1]^(N-1)*x;

% Get simulated instantaneous trilateration error
errormeasx = x_meas(1,N+1:length(td))-XX(1,N+1:length(td)); 
errormeasy = x_meas(4,N+1:length(td))-XX(4,N+1:length(td)); 

error_matricsmeasx = [mean(errormeasx), var(errormeasx), sqrt(var(errormeasx))]
error_matricsmeasy = [mean(errormeasy), var(errormeasy), sqrt(var(errormeasy))]
fprintf('-------------------------------')
% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(td)); erroresty = x(4,:)-XX(4,N+1:length(td)); 

error_matricsestxR = [mean(errorestx), var(errorestx), sqrt(var(errorestx))]
error_matricsestyR = [mean(erroresty), var(erroresty), sqrt(var(erroresty))]
fprintf('-------------------------------')

% Plot error
figure
subplot(2,1,1),
plot(td(N+1:length(td)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('simulated instantaneous trilateration error x'),
legend('simulated instantaneous trilateration error x')
subplot(2,1,2),
plot(td(N+1:length(td)), errorestx), xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),
legend('estimated error x')

figure
subplot(2,1,1),
plot(td(N+1:length(td)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('simulated instantaneous trilateration error x'),
legend('simulated instantaneous trilateration error y')
subplot(2,1,2),
plot(td(N+1:length(td)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),
legend('estimated error y')

figure
subplot(4,1,1),
plot(td(N+1:length(td)),x_meas(1,N+1:length(td)),'LineWidth',2);hold on;
plot(td(N+1:length(td)),x(1,:),'LineWidth',2);hold on;%estimate
plot(td(N+1:length(td)),XX(1,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)), 1.5*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$x[m]$','interpreter','latex');
legend('simulated instantaneous trilateration','estimated','ground truth','initialization');title('Position x')

subplot(4,1,2),
plot(td(N+1:length(td)),x(2,:),'LineWidth',2);hold on;%estimate
plot(td(N+1:length(td)),XX(2,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),-2*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$y[m]$','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity x');xlabel('t [s]')

subplot(4,1,3),
plot(td(N+1:length(td)),x_meas(4,N+1:length(td)),'LineWidth',2); hold on; 
plot(td(N+1:length(td)),x(4,:),'LineWidth',2);hold on;%estimate
plot(td(N+1:length(td)),XX(4,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),1*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$y[m]$','interpreter','latex');
legend('simulated instantaneous trilateration','estimated','ground truth','initialization');title('Position y');xlabel('t [s]')

subplot(4,1,4),
plot(td(N+1:length(td)),x(5,:),'LineWidth',2);hold on;%estimate
plot(td(N+1:length(td)),XX(5,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),-2*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$y[m]$','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity y');xlabel('t [s]')


N = 1;
% =========================  Estimation Based on R Rdot: N = 5 ================================

% Simulated simulated instantaneous trilateration y
y = NaN(6*N,length(td)-N);

% Get yN
for l = 1:1:length(td)-N
y(:,l) = getyNPV(r_sim1,rdot_sim1,r_sim2,rdot_sim2,r_sim3,rdot_sim3,l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0]*ones(1,length(td)-N);e = NaN(1,length(td)-N);

for m = 1:1:length(td)-N %time 1 t
    
for iter = 1:1:iternum
 
E = getEPV(y(:,m),x(:,m),N);
orig = norm(E);

JacoN = getJacoN(y(:,m),x(:,m),N,2);
add = pinv(JacoN)*E;

temp = x(:,m) - time_step*add;%x - 0.5.*pinv(JacoN)*E

E_temp = getEPV(y(:,m),temp,N);
after = norm(E_temp);

if norm(E_temp) < norm(E)
    x(:,m) = temp;
else
    e(m) = norm(E_temp);
    break
end

end
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0;0,1,T,0,0,0;0,0,1,0,0,0;0,0,0,1,T,0;0,0,0,0,1,T;0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(td)); erroresty = x(4,:)-XX(4,N+1:length(td)); 

error_matricsestxRRdot = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestyRRdot = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]
fprintf('-------------------------------')

% Plot error

figure
subplot(2,1,1),plot(td(N+1:length(td)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),legend('simulated instantaneous trilateration error x'),title('simulated instantaneous trilateration error x'),
subplot(2,1,2),plot(td(N+1:length(td)), errorestx), xlabel('time [s]'),ylabel('error [m]'),legend('estimated error x'),title('estimated error x'),

figure
subplot(2,1,1),plot(td(N+1:length(td)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),legend('simulated instantaneous trilateration error y'),title('simulated instantaneous trilateration error x'),
subplot(2,1,2),plot(td(N+1:length(td)), erroresty), xlabel('time [s]'),ylabel('error [m]'),legend('estimated error y'),title('estimated error x'),

figure
subplot(4,1,1),
plot(td(N+1:length(td)),x_meas(1,N+1:length(td)),'LineWidth',2);hold on;
plot(td(N+1:length(td)),x(1,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(1,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),1.5*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$x[m]$','interpreter','latex');
legend('simulated instantaneous trilateration','estimated','ground truth','initialization');title('Position x')

subplot(4,1,2),
plot(td(N+1:length(td)),x(2,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(2,N+1:length(td)),'LineWidth',2);
ylabel('$x[m/s]$','interpreter','latex');hold on,
plot(td(N+1:length(td)),-2*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
legend('estimated','ground truth','initialization');title('Velocity x')

subplot(4,1,3),
plot(td(N+1:length(td)),x_meas(4,N+1:length(td)),'LineWidth',2); hold on; 
plot(td(N+1:length(td)),x(4,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(4,N+1:length(td)),'LineWidth',2);
ylabel('$y[m]$','interpreter','latex');hold on,
plot(td(N+1:length(td)),ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
legend('simulated instantaneous trilateration','estimated','ground truth','initialization');title('Position y')

subplot(4,1,4),
plot(td(N+1:length(td)),x(5,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(5,N+1:length(td)),'LineWidth',2);
ylabel('$y[m/s]$','interpreter','latex');hold on,
plot(td(N+1:length(td)),-2*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
legend('estimated','ground truth','initialization');title('Velocity y');xlabel('td [s]')

N = 1;
% =============================  Estimation Based on R Acc ================================

% Simulated simulated instantaneous trilateration y
y = NaN(5*N,length(td)-N);

% Get yN
for l = 1:1:length(td)-N
y(:,l) = getyNPA(r_sim1,r_sim2,r_sim3,A(1,:),A(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0]*ones(1,length(td)-N);e = NaN(1,length(td)-N);

for m = 1:1:length(td)-N %time 1 t
    
for iter = 1:1:iternum
 
E = getEPA(y(:,m),x(:,m),N);
orig = norm(E);

JacoN = getJacoN(y(:,m),x(:,m),N,3);
add = pinv(JacoN)*E;

temp = x(:,m) - time_step*add;%x - 0.5.*pinv(JacoN)*E

E_temp = getEPA(y(:,m),temp,N);
after = norm(E_temp);

if norm(E_temp) < norm(E)
    x(:,m) = temp;
else
    e(m) = norm(E_temp);
    break
end

end
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0;0,1,T,0,0,0;0,0,1,0,0,0;0,0,0,1,T,0;0,0,0,0,1,T;0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(td)); erroresty = x(4,:)-XX(4,N+1:length(td)); 

error_matricsestxRAcc = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestyRAcc = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]
fprintf('-------------------------------')

% Plot error
figure
subplot(2,1,1),plot(td(N+1:length(td)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('simulated instantaneous trilateration error x'),legend('simulated instantaneous trilateration error x')
subplot(2,1,2),plot(td(N+1:length(td)), errorestx), xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),legend('estimated error x')

figure
subplot(2,1,1),plot(td(N+1:length(td)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('simulated instantaneous trilateration error y'),legend('simulated instantaneous trilateration error y')
subplot(2,1,2),plot(td(N+1:length(td)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error y'),legend('estimated error y')

figure
subplot(6,1,1),
plot(td(N+1:length(td)),x_meas(1,N+1:length(td)),'LineWidth',2);hold on;
plot(td(N+1:length(td)),x(1,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(1,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),1.5*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$x[m]$','interpreter','latex');
legend('simulated instantaneous trilateration','estimated','ground truth','initialization');title('Position x')

subplot(6,1,2),
plot(td(N+1:length(td)),x(2,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(2,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),-2*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot x$[m/s]','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity x')

subplot(6,1,3),
plot(td(N+1:length(td)),x(3,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),A(1,N+1:length(td)),'LineWidth',2);hold on,%XX(3,N+1:length(td))
plot(td(N+1:length(td)),0*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('estimated','meas','initialization');title('Acceleration x')

subplot(6,1,4),
plot(td(N+1:length(td)),x_meas(4,N+1:length(td)),'LineWidth',2); hold on; 
plot(td(N+1:length(td)),x(4,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(4,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$y[m]$','interpreter','latex');
legend('simulated instantaneous trilateration','estimated','ground truth','initialization');title('Position y')

subplot(6,1,5),
plot(td(N+1:length(td)),x(5,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(5,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),-2*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot y$[m/s]','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity y');xlabel('td [s]')

subplot(6,1,6),
plot(td(N+1:length(td)),x(6,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),A(2,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),0*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);%XX(6,N+1:length(td))
ylabel('$\ddot y[m/s^2]$','interpreter','latex');
legend('estimated','meas','initialization');title('Acceleration y');xlabel('td [s]')


N = 1;
% ======================  Estimation R Rdot Acc: N = 5 ==================================

% Simulated simulated instantaneous trilateration y
y = NaN(N*8,length(td)-N);

% Get yN
for l = 1:1:length(td)-N
y(:,l) = getyNPVA(r_sim1,rdot_sim1,r_sim2,rdot_sim2,r_sim3,rdot_sim3,A(1,:),A(2,:),l,N);
end

% Get simulated instantaneous trilateration error
errormeasx = x_meas(1,N+1:length(td))-XX(1,N+1:length(td)); 
errormeasy = x_meas(4,N+1:length(td))-XX(4,N+1:length(td)); 

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0]*ones(1,length(td)-N);e = NaN(1,length(td)-N);

for m = 1:1:length(td)-N %time 1 t
    
for iter = 1:1:iternum
    
E = getEPVA(y(:,m),x(:,m),N);
orig = norm(E);

JacoN = getJacoN(y(:,m),x(:,m),N,4);
add = pinv(JacoN)*E;

temp = x(:,m) - time_step*add;%x - 0.5.*pinv(JacoN)*E

E_temp = getEPVA(y(:,m),temp,N);
after = norm(E_temp);

if norm(E_temp) < norm(E)
    x(:,m) = x(:,m) - time_step*pinv(JacoN)*E; %modi
else
    e(m) = norm(E_temp);
    break
end

end
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0;0,1,T,0,0,0;0,0,1,0,0,0;0,0,0,1,T,0,;0,0,0,0,1,T;0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(td)); erroresty = x(4,:)-XX(4,N+1:length(td)); 

error_matricsestxRRdotAcc = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestyRRdotAcc = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

% Error
figure
subplot(2,1,1),plot(td(N+1:length(td)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),legend('simulated instantaneous trilateration error x'),title('simulated instantaneous trilateration error x'),
subplot(2,1,2),plot(td(N+1:length(td)), errorestx), xlabel('time [s]'),ylabel('error [m]'),legend('estimated error x'),title('estimated error x'),

figure
subplot(2,1,1),plot(td(N+1:length(td)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),legend('simulated instantaneous trilateration error y'),title('simulated instantaneous trilateration error x'),
subplot(2,1,2),plot(td(N+1:length(td)), erroresty), xlabel('time [s]'),ylabel('error [m]'),legend('estimated error y'),title('estimated error x'),

% Estimation Results
figure
subplot(6,1,1),
plot(td(N+1:length(td)),x_meas(1,N+1:length(td)),'LineWidth',2);hold on;
plot(td(N+1:length(td)),x(1,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(1,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),1.5*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$x[m]$','interpreter','latex');
legend('simulated instantaneous trilateration','estimated','ground truth','initialization');title('Position x')

subplot(6,1,2),
plot(td(N+1:length(td)),x(2,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(2,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),-2*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot x[m/s]$','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity x')

subplot(6,1,3),
plot(td(N+1:length(td)),x(3,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(3,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),0*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('estimated','ground truth','initialization');title('Acceleration x')

subplot(6,1,4),
plot(td(N+1:length(td)),x_meas(4,N+1:length(td)),'LineWidth',2); hold on; 
plot(td(N+1:length(td)),x(4,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(4,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$y[m/s]$','interpreter','latex');
legend('simulated instantaneous trilateration','estimated','ground truth','initialization');title('Position y')

subplot(6,1,5),
plot(td(N+1:length(td)),x(5,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(5,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),-2*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot y[m/s]$','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity y');xlabel('td [s]')

subplot(6,1,6),
plot(td(N+1:length(td)),x(6,:),'LineWidth',2);hold on;
plot(td(N+1:length(td)),XX(6,N+1:length(td)),'LineWidth',2);hold on,
plot(td(N+1:length(td)),0*ones(1,length(td)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot y[m/s^2]$','interpreter','latex');
legend('estimated','ground truth','initialization');title('Acceleration y');xlabel('td [s]')