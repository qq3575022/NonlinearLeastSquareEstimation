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

fprintf('===================== R =======================')
N = xr

% Simulated instantaneous trilateration y
y = NaN(3*N,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNP(r1f,r2f,r3f,l,N);
end

save('y.mat','y')

% Initiate x, e and temproary variable
x = [1.5;-2;0;1;-2;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N
x(:,m)=lsqnonlin(@(xx)getEP(y(:,m), xx, N),[1.5;-2;0;1;-2;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0;
     0,1,T,0,0,0;
     0,0,1,0,0,0;
     0,0,0,1,T,0;
     0,0,0,0,1,T;
     0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

fprintf('---------------min error estimation----------------')

error_matricsestR(1,:) = [sqrt(mean(errorestx)^2+mean(erroresty)^2), sqrt(var(errorestx)^2+var(erroresty)^2),sqrt(var(errorestx)+var(erroresty))]
error_matricsestRx(1,:) = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
error_matricsestRy(1,:) = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]


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
legend('estimated','ground truth','initialization');title('Velocity y');xlabel('t [s]')
