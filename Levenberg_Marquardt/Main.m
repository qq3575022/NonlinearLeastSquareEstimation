clc, clear, close all
% ============================  Load instantaneous trilateration r rdot a  ==========================
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

% ---------------------------------------------
% Get instantaneous trilateration error

% errormeasx = x_meas(1,N+1:length(tf))-XX(1,N+1:length(tf)); 
% errormeasy = x_meas(4,N+1:length(tf))-XX(4,N+1:length(tf)); 
% 
% fprintf('------------- measurement error ------------------')
% error_matricsmeas  = [sqrt(mean(errormeasx)^2+mean(errormeasy)^2), sqrt(var(errormeasx)^2+var(errormeasy)^2),sqrt(var(errormeasx) + var(errormeasy))]
% error_matricsmeasx = [mean(errormeasx), var(errormeasx),sqrt(var(errormeasx))]
% error_matricsmeasy = [mean(errormeasy), var(errormeasy),sqrt(var(errormeasy))]

% Position of readers
x1 = [-0.05, 1.5]; x2 = [2, 3.0]; x3 = [2.7, 0.05];

% Set timestep
time_step    = 0.1;   iternum    = 200; 
time_stepPos = 0.1;   iternumPos = 200;  

% ++++++++++++++++++++++   Input: Length of NLS estimation  ++++++++++++++++
meanrr = NaN(6,40);  stdrr = NaN(6,40);  time = linspace(1,40,40);
meanx = NaN(6,40);   stdx = NaN(6,40);
meany = NaN(6,40);   stdy = NaN(6,40);
%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%==============================  Estimation Based on R ===============================
for N = 1:1:40
% Simulated instantaneous trilateration y
y = NaN(3*N,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNP(r1f,r2f,r3f,l,N);
end

% Initiate x, e and temproary variable
x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N
x(:,m)=lsqnonlin(@(xx)getEP(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
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

errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 

meanx(1,N) = mean(errorestx);      meany(1,N) = mean(erroresty);
stdx(1,N) = sqrt(var(errorestx));  stdy(1,N) =  sqrt(var(erroresty));

meanrr(1,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
stdrr(1,N)  = sqrt(var(errorestx)    + var(erroresty));

end

%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% =============================  Estimation Based on R Rdot ==========================
for N = 1:1:40
% Simulated instantaneous trilateration y
y = NaN(6*N,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPV(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,l,N);
end

% Initiate x, e and temproary variable
x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPV(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
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

% error_matricsestxRRdot = [mean(errorestx), var(errorestx),sqrt(var(errorestx))]
% error_matricsestyRRdot = [mean(erroresty), var(erroresty),sqrt(var(erroresty))]
% fprintf('-------------------------------')
meanx(2,N) = mean(errorestx);     meany(2,N) = mean(erroresty);
stdx(2,N) = sqrt(var(errorestx)); stdy(2,N) =  sqrt(var(erroresty));

meanrr(2,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
stdrr(2,N)  = sqrt(var(errorestx)    + var(erroresty));
end

%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ==========================  Estimation R Acc Orientation ===========================
for N = 1:1:40
y = NaN(N*6,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPAO(r1f,r2f,r3f,phi(1,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPAO(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
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

meanx(3,N) = mean(errorestx);     meany(3,N) = mean(erroresty);
stdx(3,N) = sqrt(var(errorestx)); stdy(3,N) =  sqrt(var(erroresty));

meanrr(3,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
stdrr(3,N)  = sqrt(var(errorestx)    + var(erroresty));
end

%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% =============================  Estimation Based on R Acc ===========================
for N = 1:1:40
y =   NaN(7*N,length(tf)-N); 

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPA(r1f,r2f,r3f,phi(1,:),phi(2,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
x(:,m)=lsqnonlin(@(xx)getEPA(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
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

meanx(4,N) = mean(errorestx);     meany(4,N) = mean(erroresty);
stdx(4,N) = sqrt(var(errorestx)); stdy(4,N) =  sqrt(var(erroresty));

meanrr(4,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
stdrr(4,N)  = sqrt(var(errorestx)    + var(erroresty));

end

%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ==========================  Estimation R Rdot Acc Orientation ======================
for N = 1:1:40
y = NaN(N*9,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPVAO(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPVAO(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
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

meanx(5,N) = mean(errorestx);     meany(5,N) = mean(erroresty);
stdx(5,N) = sqrt(var(errorestx)); stdy(5,N) =  sqrt(var(erroresty));

meanrr(5,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
stdrr(5,N)  = sqrt(var(errorestx)    + var(erroresty));
end


%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ==========================  Estimation R Rdot Acc ==================================
for N = 1:1:40
y = NaN(N*10,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPVA(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),phi(2,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPVA(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
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

meanx(6,N) = mean(errorestx);     meany(6,N) = mean(erroresty);
stdx(6,N) = sqrt(var(errorestx)); stdy(6,N) =  sqrt(var(erroresty));

meanrr(6,N) = sqrt(mean(errorestx)^2 + mean(erroresty)^2); 
stdrr(6,N)  = sqrt(var(errorestx)    + var(erroresty));
end

%%
errormeasx = x_meas(1,:) - XX(1,:); errormeasy = x_meas(4,:) - XX(4,:);

time = time(1:40);
figure
subplot(2,1,1), plot(time,meanrr(1,1:40),time,meanrr(2,1:40),time,meanrr(3,1:40),time,meanrr(4,1:40),time,meanrr(5,1:40),time,meanrr(6,1:40),'LineWidth',2); hold on; plot(time, ones(1,40)*sqrt(mean(errormeasx)^2+mean(errormeasy)^2),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
grid on; xlim([1,40]);xlabel('Stack Length');ylim([0,0.03])
legend('estimated r','estimated r $\dot r$','estimated r $a_x$ $a_y$ $\theta_z$','estimated r $a_x$ $a_y$ $\theta_z$ $\omega_z$','estimated r $\dot r$ $a_x$ $a_y$ $\theta_z$','estimated r $\dot r$ $a_x$ $a_y$ $\theta_z$ $\omega_z$','interpreter','latex','Location','NorthWest');title('Mean Error of Estimation');ylabel('Mean Error [m]');xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])

subplot(2,1,2), plot(time,stdrr(1,1:40),time,stdrr(2,1:40),time,stdrr(3,1:40),time,stdrr(4,1:40),time,stdrr(5,1:40),time,stdrr(6,1:40),'LineWidth',2);hold on; plot(time, ones(1,40)*sqrt(var(errormeasx) + var(errormeasy)),'Color',[0.4940 0.1840 0.5560],'LineStyle',':','LineWidth',2);
grid on; xlim([1,40]);xlabel('Stack Length');ylim([0.015,0.075]);
legend('estimated r','estimated r $\dot r$','estimated r $a_x$ $a_y$ $\theta_z$','estimated r $a_x$ $a_y$ $\theta_z$ $\omega_z$','estimated r $\dot r$ $a_x$ $a_y$ $\theta_z$','estimated r $\dot r$ $a_x$ $a_y$ $\theta_z$ $\omega_z$','instantaneous trilateration','interpreter','latex','Location','NorthWest');
title('RMS Error of Estimation');ylabel('RMS Error [m]');xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])

%%

save('meanrr.mat','meanrr')
save('stdrr.mat', 'stdrr')

save('meanx_3.mat','meanx')
save('stdx_3.mat','stdx')

save('meany_3.mat','meany')
save('stdy_3.mat','stdy')

%%
fprintf('----------- E ----------')
meanE = sqrt(mean(x_meas(1,:) - XX(1,:))^2 + mean(x_meas(4,:) - XX(4,:))^2)
varE = sqrt(var(x_meas(1,:) - XX(1,:)) + var(x_meas(4,:) - XX(4,:)))

fprintf('----------- R ----------')
meanR = meanrr(1,10)
minxR = min(meanrr(1,:))
meanR = find(meanrr(1,:)==minxR)

varR  = stdrr(1,10)
minyR = min(stdrr(1,:)) 
varR  = find(stdrr(1,:)==minyR)

fprintf('----------- R R----------')
meanRR = meanrr(2,10)
minxRR = min(meanrr(2,:))
meanRR = find(meanrr(2,:)==minxRR)

varRR  = stdrr(2,10)
minyRR = min(stdrr(2,:)) 
varRR  = find(stdrr(2,:)==minyRR)

fprintf('----------- R ACC O----------')
meanRAO = meanrr(3,10)
minxRAO = min(meanrr(3,:))
meanRAO = find(meanrr(3,:)==minxRAO)

varRAO  = stdrr(3,10)
minyRAO = min(stdrr(3,:)) 
varRAO  = find(stdrr(3,:)==minyRAO)

fprintf('----------- R ACC ----------')
meanRA = meanrr(4,10)
minxRA = min(meanrr(4,:))
meanRA = find(meanrr(4,:)==minxRA)

varRA  = stdrr(4,10)
minyRA = min(stdrr(4,:)) 
varRA  = find(stdrr(4,:)==minyRA)

fprintf('----------- R R ACC O----------')
meanRRAO = meanrr(5,10)
minxRRAO = min(meanrr(5,:))
meanRRAO = find(meanrr(5,:)==minxRRAO)

varRRAO  = stdrr(5,10)
minyRRAO = min(stdrr(5,:))  
varRRAO  = find(stdrr(5,:)==minyRRAO)

fprintf('----------- R R ACC ----------')
meanRRA = meanrr(6,10)
minxRRA = min(meanrr(6,:)) 
meanRRA = find(meanrr(6,:)==minxRRA)

varRRA  = stdrr(6,10)
minyRRA = min(stdrr(6,:))  
varRRA  = find(stdrr(6,:)==minyRRA)

