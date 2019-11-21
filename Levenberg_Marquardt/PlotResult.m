close all
load('yVectorData.mat');global T; T = mean(diff(tf));
% -----------Load Data Infromation------------
% time: tf    
% radical distances: r1f, r2f, r3f
% radical velocity:  r1dot_e, r2dot_e, r3dot_e
% orientation    psif1
% angular velocity  wf
% acceleration ax1 ay1
% ---------------------------------------------
% Ground truth coordinates of the moving tag path

load('meanrr.mat')
load('stdrr.mat')


% Set timestep
time_step    = 0.1;    iternum    = 200; 
time_stepPos = 0.1;    iternumPos = 200;  
% ---------------------------------------------
% Get acceleration along x^B y^B
axx = ax1+0.6;  axx = 0.8*axx;  ayy = ay1+0.15; ayy = 0.8*ayy;
AXY = NaN(2,length(axx));       AXY(1,:) = axx; AXY(2,:) = ayy;

% ---------------------------------------------
% Get measurement orientation and angular velocity
offset = 0;
for n = 2:1:length(psif1)
    if psif1(n) - psif1(n-1) < -2
        offset = 2*pi;
    end
    psif1(n) = psif1(n) + offset;  
end
phi = NaN(2,length(axx)); phi(1,:) = psif1 - pi; phi(2,:) = 0.85*wf;

% ---------------------------------------------
% Instantaneous Trilateration
[x_meas, A] = rtoxy(r1f, r2f, r3f, r1dot_e, r2dot_e, r3dot_e, AXY, T);

% ---------------------------------------------
% Get ground-truth
[W_A, W_V, W, XX, rr, AXY_gt, Orient] = loadgtruth(tf);

errormeasx = x_meas(1,:) - XX(1,:); errormeasy = x_meas(4,:) - XX(4,:);

%----------------------  Estimation Based on R min Error------------------------
minx = min(meanrr(1,:)); xr =find(meanrr(1,:)==minx)
miny = min(stdrr(1,:));  yr =find(stdrr(1,:)==miny)

errormeasx = x_meas(1,:) - XX(1,:); errormeasy = x_meas(4,:) - XX(4,:);

figure
subplot(2,1,1), plot(time,meanrr(1,1:40),time,meanrr(2,1:40),time,meanrr(3,1:40),time,meanrr(4,1:40),time,meanrr(5,1:40),time,meanrr(6,1:40),'LineWidth',2); hold on; plot(time, ones(1,40)*sqrt(mean(errormeasx)^2+mean(errormeasy)^2),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
grid on; xlim([1,40]);xlabel('Stack Length');ylim([0,0.03])
legend('estimated r','estimated r rdot','estimated r acc $\phi$','estimated r acc $\phi$ $\dot \phi$','estimated r rdot acc $\phi$','estimated r rdot acc $\phi$ $\dot \phi$','instantaneous trilateration','interpreter','latex');title('error mean of estimation');ylabel('error mean [m]');xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])

subplot(2,1,2), plot(time,stdrr(1,1:40),time,stdrr(2,1:40),time,stdrr(3,1:40),time,stdrr(4,1:40),time,stdrr(5,1:40),time,stdrr(6,1:40),'LineWidth',2);hold on; plot(time, ones(1,40)*sqrt(var(errormeasx) + var(errormeasy)),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
grid on; xlim([1,40]);xlabel('Stack Length');ylim([0.015,0.075]);
legend('estimated r','estimated r rdot','estimated r acc $\phi$','estimated r acc $\phi$ $\dot \phi$','estimated r rdot acc $\phi$','estimated r rdot acc $\phi$ $\dot \phi$','instantaneous trilateration','interpreter','latex');title('error standard deviation of estimation');ylabel('error sd [m]');xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])
%%
%meanrr
%stdrr
fprintf('===================== R =======================')
global N
N = xr

% Simulated instantaneous trilateration y
y = NaN(3*N,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNP(r1f,r2f,r3f,l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

global m
for m = 1:1:length(tf)-N
x(:,m)=lsqnonlin(@(xx)getEP(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('---------------min error estimation----------------')

error_matricsestR(1,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRx(1,:) = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRy(1,:) = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

% Plot error
% figure
% subplot(2,1,1),
% plot(tf(N+1:length(tf)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error x'),
% legend('instantaneous trilateration error x');
% subplot(2,1,2),
% plot(tf(N+1:length(tf)), errorestx), xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),
% legend('estimated error x')

% figure
% subplot(2,1,1),
% plot(tf(N+1:length(tf)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error y'),
% legend('instantaneous trilateration error y')
% subplot(2,1,2),
% plot(tf(N+1:length(tf)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error y'),
% legend('estimated error y')

figure
subplot(4,1,1),
plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)), 1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$x[m]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position x')

subplot(4,1,2),
plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$\dot x[m]$','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity x');xlabel('t [s]')

subplot(4,1,3),
plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),1*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$y[m]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position y');xlabel('t [s]')

subplot(4,1,4),
plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;%estimate
plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
ylabel('$\dot y[m]$','interpreter','latex');
legend({'estimated','ground truth','initialization'},'Location','NorthEast');title('Velocity y');xlabel('t [s]')

N = yr
% Simulated instantaneous trilateration y
y = NaN(3*N,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNP(r1f,r2f,r3f,l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N
x(:,m)=lsqnonlin(@(xx)getEP(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('---------------min variance estimation----------------')

error_matricsestR(2,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRx(2,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRy(2,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]
% Plot error
% figure
% subplot(2,1,1),
% plot(tf(N+1:length(tf)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error x'),
% legend('instantaneous trilateration error x');
% subplot(2,1,2),
% plot(tf(N+1:length(tf)), errorestx), xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),
% legend('estimated error x')

% figure
% subplot(2,1,1),
% plot(tf(N+1:length(tf)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error y'),
% legend('instantaneous trilateration error y')
% subplot(2,1,2),
% plot(tf(N+1:length(tf)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error y'),
% legend('estimated error y')

% figure
% subplot(4,1,1),
% plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;%estimate
% plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)), 1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
% ylabel('$x[m]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position x')
% 
% subplot(4,1,2),
% plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;%estimate
% plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
% ylabel('$\dot x[m]$','interpreter','latex');
% legend('estimated','ground truth','initialization');title('Velocity x');xlabel('t [s]')
% 
% subplot(4,1,3),
% plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
% plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;%estimate
% plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),1*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
% ylabel('$y[m]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position y');xlabel('t [s]')
% 
% subplot(4,1,4),
% plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;%estimate
% plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2)
% ylabel('$\dot y[m]$','interpreter','latex');
% legend('estimated','ground truth','initialization');title('Velocity y');xlabel('t [s]')

N = 10
y = NaN(3*N,length(tf)-N);
% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNP(r1f,r2f,r3f,l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N
x(:,m)=lsqnonlin(@(xx)getEP(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('---------------length = 10 estimation----------------')

error_matricsestR(3,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRx(3,:) = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRy(3,:) = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

%%
fprintf('======================= R Rdot ==========================')
% =============================  Estimation Based on R Rdot ================================
minx = min(meanrr(2,:)); xr =find(meanrr(2,:)==minx)
miny = min(stdrr(2,:));  yr =find(stdrr(2,:)==miny)

N = xr
% Simulated instantaneous trilateration y
y = NaN(6*N,length(tf)-N);
% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPV(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
x(:,m)=lsqnonlin(@(xx)getEPV(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('---------------min error estimation----------------')
error_matricsestxRR(1,:)  = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestyRRx(1,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestyRRy(1,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]
%%
% Plot error
figure
subplot(4,1,1),
plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$x[m]$','interpreter','latex');
legend({'Instantaneous trilateration','Estimated','Ground truth','Initialization'},'Location','SouthEast');title('Position $x$','interpreter','latex');xlabel('t[s]')

subplot(4,1,2),
plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',2);
ylabel('$\dot x[m/s]$','interpreter','latex');hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
legend('Estimated','Ground truth','Initialization');title('Velocity $\dot x$','interpreter','latex');xlabel('t[s]')

subplot(4,1,3),
plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',2);
ylabel('$y[m]$','interpreter','latex');hold on,
plot(tf(N+1:length(tf)),ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
legend({'Instantaneous trilateration','Estimated','Ground truth','Initialization'},'Location','SouthEast');title('Position $y$','interpreter','latex');xlabel('t[s]')

subplot(4,1,4),
plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);
ylabel('$\dot y[m/s]$','interpreter','latex');hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
legend('Estimated','Ground truth','Initialization');title('Velocity $\dot y$','interpreter','latex');xlabel('t[s]')

%%
N = yr
% Simulated instantaneous trilateration y
y = NaN(6*N,length(tf)-N);
% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPV(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPV(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('------------- mean error variance estimation ------------------')
error_matricsestxRR(2,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestyRRx(2,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestyRRy(2,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

% Plot error
% figure
% subplot(4,1,1),
% plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$x[m]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position x')
% 
% subplot(4,1,2),
% plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',2);
% ylabel('$\dot x[m/s]$','interpreter','latex');hold on,
% plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% legend('estimated','ground truth','initialization');title('Velocity x')
% 
% subplot(4,1,3),
% plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
% plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',2);
% ylabel('$y[m]$','interpreter','latex');hold on,
% plot(tf(N+1:length(tf)),ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position y')
% 
% subplot(4,1,4),
% plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);
% ylabel('$\dot y[m/s]$','interpreter','latex');hold on,
% plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% legend('estimated','ground truth','initialization');title('Velocity y');xlabel('tf [s]')


N = 10;
y = NaN(6*N,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPV(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
x(:,m)=lsqnonlin(@(xx)getEPV(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('----------- N = 10 estimation --------------------')
error_matricsestRR(3,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRRx(3,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRRy(3,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

%%
% ==========================  Estimation R Acc Orientation ==================================
fprintf('======================= R Acc Orientation  ==========================')
minx = min(meanrr(3,:)); xr =find(meanrr(3,:)==minx)
miny = min(stdrr(3,:));  yr =find(stdrr(3,:)==miny)

N = xr%min(x,y)

y = NaN(N*6,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPAO(r1f,r2f,r3f,phi(1,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPAO(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('------------- mean error estimation ------------------')
error_matricsestRAO(1,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRAOx(1,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRAOy(1,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

% Error
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error x'),legend('instantaneous trilateration error x')
% subplot(2,1,2),plot(tf(N+1:length(tf)), errorestx),xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),legend('estimated error x')
% 
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error y'),legend('instantaneous trilateration error y')
% subplot(2,1,2),plot(tf(N+1:length(tf)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error y'),legend('estimated error y')

% Estimation Results
figure
subplot(8,1,1),
plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$x[m]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position x')

subplot(8,1,2),
plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot x$[m/s]','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity x')

subplot(8,1,3),
plot(tf(N+1:length(tf)),A(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(3,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(3,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('instantaneous trilateration','estimated','meas','initialization');title('Acceleration x')

subplot(8,1,4),
plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$y[m]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position y')

subplot(8,1,5),
plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot y$[m/s]','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity y');xlabel('tf [s]')

subplot(8,1,6),
plot(tf(N+1:length(tf)),A(2,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(6,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(6,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);%XX(6,N+1:length(tf))
ylabel('$\ddot y[m/s^2]$','interpreter','latex');
legend('instantaneous trilateration','estimated','meas','initialization');title('Acceleration y');xlabel('tf [s]')

subplot(8,1,7),
plot(tf(N+1:length(tf)),phi(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),x(7,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),W(N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('meas','estimated','ground truth','initialization');title('Orientation')

subplot(8,1,8),
plot(tf(N+1:length(tf)),phi(2,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),x(8,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),W_V(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('meas','estimated','ground truth','initialization');title('Angular Velocity')

N = yr%min(x,y)

y = NaN(N*6,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPAO(r1f,r2f,r3f,phi(1,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPAO(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('------------- mean error variance estimation ------------------')
error_matricsestRAO(2,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRAOx(2,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRAOy(2,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

% Error
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error x'),legend('instantaneous trilateration error x')
% subplot(2,1,2),plot(tf(N+1:length(tf)), errorestx),xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),legend('estimated error x')
% 
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error y'),legend('instantaneous trilateration error y')
% subplot(2,1,2),plot(tf(N+1:length(tf)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error y'),legend('estimated error y')

% Estimation Results
% figure
% subplot(8,1,1),
% plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$x[m]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position x')
% 
% subplot(8,1,2),
% plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\dot x$[m/s]','interpreter','latex');
% legend('estimated','ground truth','initialization');title('Velocity x')
% 
% subplot(8,1,3),
% plot(tf(N+1:length(tf)),A(1,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(3,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(3,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\ddot x[m/s^2]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','meas','initialization');title('Acceleration x')
% 
% subplot(8,1,4),
% plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
% plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$y[m]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position y')
% 
% subplot(8,1,5),
% plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\dot y$[m/s]','interpreter','latex');
% legend('estimated','ground truth','initialization');title('Velocity y');xlabel('tf [s]')
% 
% subplot(8,1,6),
% plot(tf(N+1:length(tf)),A(2,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(6,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(6,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);%XX(6,N+1:length(tf))
% ylabel('$\ddot y[m/s^2]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','meas','initialization');title('Acceleration y');xlabel('tf [s]')
% 
% subplot(8,1,7),
% plot(tf(N+1:length(tf)),phi(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),x(7,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),W(N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\ddot x[m/s^2]$','interpreter','latex');
% legend('meas','estimated','ground truth','initialization');title('Orientation')
% 
% subplot(8,1,8),
% plot(tf(N+1:length(tf)),phi(2,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),x(8,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),W_V(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\ddot x[m/s^2]$','interpreter','latex');
% legend('meas','estimated','ground truth','initialization');title('Angular Velocity')

N = 10
y = NaN(N*6,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPAO(r1f,r2f,r3f,phi(1,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPAO(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('------------- N = 10 estimation ------------------')
error_matricsestRAO(3,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRAOx(3,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRAOy(3,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

%%
% =============================  Estimation Based on R Acc ================================
fprintf('======================= R Acc ==========================')
minx = min(meanrr(4,:)); xr =find(meanrr(4,:)==minx)
miny = min(stdrr(4,:));  yr =find(stdrr(4,:)==miny)


N = xr%min(x,y)
y =   NaN(7*N,length(tf)-N); 
% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPA(r1f,r2f,r3f,phi(1,:),phi(2,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
x(:,m)=lsqnonlin(@(xx)getEPA(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('------------- mean error estimation ------------------')
error_matricsestRA(1,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRAx(1,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRAy(1,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error x'),legend('instantaneous trilateration error x')
% subplot(2,1,2),plot(tf(N+1:length(tf)), errorestx), xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),legend('estimated error x')
% 
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error y'),legend('instantaneous trilateration error y')
% subplot(2,1,2),plot(tf(N+1:length(tf)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error y'),legend('estimated error y')

figure
subplot(8,1,1),
plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',1.5);hold on,
plot(tf(N+1:length(tf)),1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$x[m]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position x')

subplot(8,1,2),
plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',1.5);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot x$[m/s]','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity x')

subplot(8,1,3),
plot(tf(N+1:length(tf)),A(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(3,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(3,N+1:length(tf)),'LineWidth',1.5);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Acceleration x')

subplot(8,1,4),
plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',1.5);hold on,
plot(tf(N+1:length(tf)),ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$y[m]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position y')

subplot(8,1,5),
plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot y$[m/s]','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity y');xlabel('tf [s]')

subplot(8,1,6),
plot(tf(N+1:length(tf)),A(2,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(6,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(6,N+1:length(tf)),'LineWidth',1.5);hold on,
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);%XX(6,N+1:length(tf))
ylabel('$\ddot y[m/s^2]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Acceleration y');xlabel('tf [s]')

subplot(8,1,7),
plot(tf(N+1:length(tf)),phi(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),x(7,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),Orient(1,N+1:length(tf)),'LineWidth',1.5);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('meas','estimated','ground truth','initialization');title('Orientation')

subplot(8,1,8),
plot(tf(N+1:length(tf)),phi(2,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),x(8,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),W_V(N+1:length(tf)),'LineWidth',1.5);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('meas','estimated','ground truth','initialization');title('Angular Velocity')

N = yr%min(x,y)
y =   NaN(7*N,length(tf)-N); 
% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPA(r1f,r2f,r3f,phi(1,:),phi(2,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
x(:,m)=lsqnonlin(@(xx)getEPA(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('------------- mean error variance estimation ------------------')
error_matricsestRA(2,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRAx(2,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRAy(2,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error x'),legend('instantaneous trilateration error x')
% subplot(2,1,2),plot(tf(N+1:length(tf)), errorestx), xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),legend('estimated error x')
% 
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error y'),legend('instantaneous trilateration error y')
% subplot(2,1,2),plot(tf(N+1:length(tf)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error y'),legend('estimated error y')

% figure
% subplot(8,1,1),
% plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',1.5);hold on,
% plot(tf(N+1:length(tf)),1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$x[m]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position x')
% 
% subplot(8,1,2),
% plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',1.5);hold on,
% plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\dot x$[m/s]','interpreter','latex');
% legend('estimated','ground truth','initialization');title('Velocity x')
% 
% subplot(8,1,3),
% plot(tf(N+1:length(tf)),A(1,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(3,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(3,N+1:length(tf)),'LineWidth',1.5);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\ddot x[m/s^2]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Acceleration x')
% 
% subplot(8,1,4),
% plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
% plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',1.5);hold on,
% plot(tf(N+1:length(tf)),ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$y[m]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position y')
% 
% subplot(8,1,5),
% plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\dot y$[m/s]','interpreter','latex');
% legend('estimated','ground truth','initialization');title('Velocity y');xlabel('tf [s]')
% 
% subplot(8,1,6),
% plot(tf(N+1:length(tf)),A(2,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(6,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(6,N+1:length(tf)),'LineWidth',1.5);hold on,
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);%XX(6,N+1:length(tf))
% ylabel('$\ddot y[m/s^2]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Acceleration y');xlabel('tf [s]')
% 
% subplot(8,1,7),
% plot(tf(N+1:length(tf)),phi(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),x(7,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),Orient(1,N+1:length(tf)),'LineWidth',1.5);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\ddot x[m/s^2]$','interpreter','latex');
% legend('meas','estimated','ground truth','initialization');title('Orientation')
% 
% subplot(8,1,8),
% plot(tf(N+1:length(tf)),phi(2,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),x(8,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),W_V(N+1:length(tf)),'LineWidth',1.5);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\ddot x[m/s^2]$','interpreter','latex');
% legend('meas','estimated','ground truth','initialization');title('Angular Velocity')


N = 10
y =   NaN(7*N,length(tf)-N); 

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPA(r1f,r2f,r3f,phi(1,:),phi(2,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
x(:,m)=lsqnonlin(@(xx)getEPA(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('------------ N = 10 estimation -------------------')
error_matricsestRA(3,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRAx(3,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRAy(3,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

%%
% ==========================  Estimation R Rdot Acc Orientation ==================================
fprintf('======================= R Rdot Acc Orientation ==========================')

minx = min(meanrr(5,:)); xr =find(meanrr(5,:)==minx)
miny = min(stdrr(5,:));  yr =find(stdrr(5,:)==miny)

N = xr

y = NaN(N*9,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPVAO(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPVAO(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('---------- min error estimation ---------------------')
error_matricsestRRAO(1,:)  = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRRAOx(1,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRRAOy(1,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

% Error
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error x'),legend('instantaneous trilateration error x')
% subplot(2,1,2),plot(tf(N+1:length(tf)), errorestx),xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),legend('estimated error x')
% 
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error y'),legend('instantaneous trilateration error y')
% subplot(2,1,2),plot(tf(N+1:length(tf)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error y'),legend('estimated error y')

% Estimation Results
figure
subplot(8,1,1),
plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$x[m]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position x')

subplot(8,1,2),
plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot x$[m/s]','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity x')

subplot(8,1,3),
plot(tf(N+1:length(tf)),A(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(3,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(3,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Acceleration x')

subplot(8,1,4),
plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$y[m]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position y')

subplot(8,1,5),
plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot y$[m/s]','interpreter','latex');
legend('estimated','ground truth','initialization');title('Velocity y');xlabel('tf [s]')

subplot(8,1,6),
plot(tf(N+1:length(tf)),A(2,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(6,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(6,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);%XX(6,N+1:length(tf))
ylabel('$\ddot y[m/s^2]$','interpreter','latex');
legend('instantaneous trilateration','estimated','ground truth','initialization');title('Acceleration y');xlabel('tf [s]')

subplot(8,1,7),
plot(tf(N+1:length(tf)),phi(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),x(7,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),W(N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('meas','estimated','ground truth','initialization');title('Orientation')

subplot(8,1,8),
plot(tf(N+1:length(tf)),phi(2,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),x(8,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),W_V(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
legend('meas','estimated','ground truth','initialization');title('Angular Velocity')

N = yr

y = NaN(N*9,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPVAO(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPVAO(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('---------- min error variance estimation ---------------------')
error_matricsestRRAO(2,:)  = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRRAOx(2,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRRAOy(2,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

% Error
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error x'),legend('instantaneous trilateration error x')
% subplot(2,1,2),plot(tf(N+1:length(tf)), errorestx),xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),legend('estimated error x')
% 
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error y'),legend('instantaneous trilateration error y')
% subplot(2,1,2),plot(tf(N+1:length(tf)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error y'),legend('estimated error y')

% Estimation Results
% figure
% subplot(8,1,1),
% plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$x[m]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position x')
% 
% subplot(8,1,2),
% plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\dot x$[m/s]','interpreter','latex');
% legend('estimated','ground truth','initialization');title('Velocity x')
% 
% subplot(8,1,3),
% plot(tf(N+1:length(tf)),A(1,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(3,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(3,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\ddot x[m/s^2]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Acceleration x')
% 
% subplot(8,1,4),
% plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
% plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$y[m]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Position y')
% 
% subplot(8,1,5),
% plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\dot y$[m/s]','interpreter','latex');
% legend('estimated','ground truth','initialization');title('Velocity y');xlabel('tf [s]')
% 
% subplot(8,1,6),
% plot(tf(N+1:length(tf)),A(2,N+1:length(tf)),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),x(6,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),XX(6,N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);%XX(6,N+1:length(tf))
% ylabel('$\ddot y[m/s^2]$','interpreter','latex');
% legend('instantaneous trilateration','estimated','ground truth','initialization');title('Acceleration y');xlabel('tf [s]')
% 
% subplot(8,1,7),
% plot(tf(N+1:length(tf)),phi(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),x(7,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),W(N+1:length(tf)),'LineWidth',2);hold on,
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\ddot x[m/s^2]$','interpreter','latex');
% legend('meas','estimated','ground truth','initialization');title('Orientation')
% 
% subplot(8,1,8),
% plot(tf(N+1:length(tf)),phi(2,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),x(8,:),'LineWidth',2);hold on;
% plot(tf(N+1:length(tf)),W_V(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
% plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
% ylabel('$\ddot x[m/s^2]$','interpreter','latex');
% legend('meas','estimated','ground truth','initialization');title('Angular Velocity')
% 

N = 10
y = NaN(N*9,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPVAO(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPVAO(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

errormeasx = x_meas(1,N+1:length(tf))-XX(1,N+1:length(tf)); 
errormeasy = x_meas(4,N+1:length(tf))-XX(4,N+1:length(tf)); 

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('------------ N = 10 estimation -------------------')
error_matricsestRRAO(3,:)  = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRRAOx(3,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRRAOy(3,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

%%
% ==========================  Estimation R Rdot Acc ==================================
fprintf('======================= R Rdot Acc ==========================')
minx = min(meanrr(6,:)); xr =find(meanrr(6,:)==minx)
miny = min(stdrr(6,:));  yr =find(stdrr(6,:)==miny)

N = xr

y = NaN(N*10,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPVA(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),phi(2,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPVA(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('----------- min error estimation --------------------')
error_matricsestRRA(1,:)  = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRRAx(1,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRRAy(1,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

% Error
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasx),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error x'),legend('instantaneous trilateration error x')
% subplot(2,1,2),plot(tf(N+1:length(tf)), errorestx), xlabel('time [s]'),ylabel('error [m]'),title('estimated error x'),legend('estimated error x')
% 
% figure
% subplot(2,1,1),plot(tf(N+1:length(tf)), errormeasy),xlabel('time [s]'),ylabel('error [m]'),title('instantaneous trilateration error y'),legend('instantaneous trilateration error y')
% subplot(2,1,2),plot(tf(N+1:length(tf)), erroresty), xlabel('time [s]'),ylabel('error [m]'),title('estimated error y'),legend('estimated error y')
%%
% Estimation Results
figure
subplot(8,1,1),
plot(tf(N+1:length(tf)),x_meas(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(1,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(1,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),1.5*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$x[m]$','interpreter','latex');
xlim([0, 2.5]);xlabel('t [s]')
legend('Instantaneous trilateration','Estimated','Ground truth','Initialization');title('Position $x$','interpreter','latex')

subplot(8,1,2),
plot(tf(N+1:length(tf)),x(2,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(2,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot x$[m/s]','interpreter','latex');
xlim([0, 2.3]);;xlabel('t [s]')
legend('Estimated','Ground truth','Initialization');title('Velocity $\dot x$','interpreter','latex')

subplot(8,1,3),
plot(tf(N+1:length(tf)),A(1,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(3,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(3,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\ddot x[m/s^2]$','interpreter','latex');
xlim([0, 2.5]);;xlabel('t [s]')
legend('Instantaneous trilateration','Estimated','Ground truth','Initialization');title('Acceleration $\ddot x$','interpreter','latex')

subplot(8,1,4),
plot(tf(N+1:length(tf)),x_meas(4,N+1:length(tf)),'LineWidth',2); hold on; 
plot(tf(N+1:length(tf)),x(4,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(4,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$y[m]$','interpreter','latex');
xlim([0, 2.5]);
legend('Instantaneous trilateration','Estimated','Ground truth','Initialization');title('Position $y$','interpreter','latex');xlabel('t [s]')

subplot(8,1,5),
plot(tf(N+1:length(tf)),x(5,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(5,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),-2*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot y$[m/s]','interpreter','latex');
xlim([0, 2.3]);
legend('Estimated','Ground truth','Initialization');title('Velocity $\dot y$','interpreter','latex');xlabel('t [s]')

subplot(8,1,6),
plot(tf(N+1:length(tf)),A(2,N+1:length(tf)),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),x(6,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),XX(6,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);%XX(6,N+1:length(tf))
ylabel('$\ddot y[m/s^2]$','interpreter','latex');
xlim([0, 2.5]);
legend('Instantaneous trilateration','Estimated','Ground truth','Initialization');title('Acceleration $\ddot y$','interpreter','latex');xlabel('t [s]')

subplot(8,1,7),
plot(tf(N+1:length(tf)),phi(1,N+1:length(tf)),'LineWidth',2);hold on,
plot(tf(N+1:length(tf)),x(7,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),W(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);xlabel('t [s]')
ylabel('$\psi[rad]$','interpreter','latex');
xlim([0, 2.2]);
legend('Measurement \theta_z','Estimated','Ground truth','Initialization');title('Orientation $\psi$','interpreter','latex')

subplot(8,1,8),
plot(tf(N+1:length(tf)),phi(2,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),x(8,:),'LineWidth',2);hold on;
plot(tf(N+1:length(tf)),W_V(1,N+1:length(tf)),'LineWidth',2);hold on,%XX(3,N+1:length(tf))
plot(tf(N+1:length(tf)),0*ones(1,length(tf)-N),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
ylabel('$\dot \psi[rad/s]$','interpreter','latex');
xlim([0, 2.2]);xlabel('t [s]')
legend('Measurement \omega_z','Estimated','Ground truth','Initialization');title('Angular Velocity $\dot \psi$','interpreter','latex')

%%
N = yr
y = NaN(N*10,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPVA(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),phi(2,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPVA(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('---------- min error variance estimation ---------------------')
error_matricsestRRA(2,:)  = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRRAx(2,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRRAy(2,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]

N = 10
y = NaN(N*10,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPVA(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),phi(2,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPVA(y(:,m), xx, N),[1.5;-2;0;1;-2;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('------------ N = 10 estimation -------------------')
error_matricsestRRA(3,:)  = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRRAx(3,:)  = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRRAy(3,:)  = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]