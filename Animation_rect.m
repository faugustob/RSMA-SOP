clc; clear; close all;

choice = menu('Select Animation', ...
    '1) rect * rect -> triangle', ...
    '2) Build periodic triangular message', ...
    '3) Spectrum of periodic triangular message', ...
    '4) AM modulation (true time-building animation)', ...
    '5) Modulation index effect (under → over)', ...
    '6) AM spectrum shifting', ...
    '7) DSB-SC phase comparison', ...
    '8) Continuous phase change in DSB-SC', ...
    '9) AM vs DSB-SC comparison');

switch choice

%% =========================================================
case 1   % Convolution -> triangle

dt = 0.002;
t = -2:dt:2;
rect1 = double(abs(t)<=0.5);

figure('Color','w')

subplot(2,1,1)
h_rect2 = plot(t,zeros(size(t)),'r','LineWidth',2); hold on
plot(t,rect1,'b','LineWidth',2)
h_overlap = area(t,zeros(size(t)),'FaceAlpha',0.3,'EdgeColor','none');
ylim([-0.1 1.2]); grid on
title('Sliding rectangles')

subplot(2,1,2)
h_conv = animatedline('LineWidth',2);
ylim([0 1.2]); grid on
title('Triangle forming from convolution')

for k=1:length(t)
    shift=t(k);
    rect2=double(abs(t-shift)<=0.5);
    overlap=rect1.*rect2;
    conv_val=sum(overlap)*dt;

    set(h_rect2,'YData',rect2)
    set(h_overlap,'YData',overlap)
    addpoints(h_conv,t(k),conv_val)
    drawnow
end


%% =========================================================
case 2   % Build periodic triangular message

dt=0.0005;
t=-0.4:dt:0.4;
T0=0.12; alpha=0.04;

figure('Color','w')
h=animatedline('LineWidth',2);
ylim([0 1.2]); grid on
title('Building periodic triangular message m(t)')

m=zeros(size(t));

for n=-5:5
    tau=t-n*T0;
    tri=max(1-abs(tau)/alpha,0);
    m=m+tri;

    clearpoints(h)
    for k=1:length(t)
        addpoints(h,t(k),m(k))
    end
    drawnow
    pause(0.4)
end


%% =========================================================
case 3   % Spectrum build-up

alpha=40e-3;
T0=120e-3;
f0=1/T0;
nmax=10;

figure('Color','w')
f=-nmax*f0:f0/100:nmax*f0;
env=alpha*f0*sinc(alpha*f).^2;

plot(f,env,'k--','LineWidth',1.5); hold on
grid on
title('Spectrum build-up')
ylim([0 max(env)*1.2])

for n=-nmax:nmax
    cn=alpha*f0*sinc(n*alpha*f0)^2;
    stem(n*f0,cn,'b','LineWidth',2)
    drawnow
    pause(0.3)
end


%% =========================================================
case 4   % AM time-building animation (correct)

dt=0.0005;
t=0:dt:0.4;
fc=100; T0=0.12; alpha=0.04;

m=zeros(size(t));
for n=-5:5
    tau=t-n*T0;
    m=m+max(1-abs(tau)/alpha,0);
end

figure('Color','w')

subplot(3,1,1)
h_m=animatedline('LineWidth',2);
ylim([0 1.2]); grid on
title('Message m(t)')

subplot(3,1,2)
h_c=animatedline('LineWidth',2);
ylim([-1.2 1.2]); grid on
title('Carrier')

subplot(3,1,3)
h_s=animatedline('LineWidth',2); hold on
plot(t,1+m,'k--')
plot(t,-(1+m),'k--')
ylim([-2 2]); grid on
title('AM signal')

for k=1:length(t)
    carrier=cos(2*pi*fc*t(k));
    s=(1+m(k))*carrier;

    addpoints(h_m,t(k),m(k))
    addpoints(h_c,t(k),carrier)
    addpoints(h_s,t(k),s)

    drawnow
end


%% =========================================================
case 5   % Modulation index sweep

dt=0.0005;
t=0:dt:0.3;
fc=100;
m=cos(2*pi*5*t);

figure('Color','w')
ylim([-2 2]); grid on
title('Modulation index sweep')

h=animatedline('LineWidth',2);

for mu=0:0.2:1.5
    clearpoints(h)
    for k=1:length(t)
        s=(1+mu*m(k))*cos(2*pi*fc*t(k));
        addpoints(h,t(k),s)
        drawnow
    end
    pause(0.5)
end


%% =========================================================
case 6   % AM spectrum shifting

alpha=40e-3;
T0=120e-3;
f0=1/T0;
fc=100;
nmax=5;

n=-nmax:nmax;
f_disc=n*f0;
M=alpha*f0*sinc(alpha*f_disc).^2;

figure('Color','w')
grid on
title('AM spectrum shifting')

for shift=0:5:fc
    cla
    stem([-fc fc],[0.5 0.5],'r','LineWidth',2); hold on
    stem(f_disc+shift,0.5*M,'b','LineWidth',2)
    stem(f_disc-shift,0.5*M,'b','LineWidth',2)
    xlim([-fc-nmax*f0-10 fc+nmax*f0+10])
    ylim([0 max(M)*0.7])
    grid on
    drawnow
end


%% =========================================================
case 7

clc; close all;

dt = 0.0005;
t = -0.1:dt:0.1;

fm = 5;
fc = 50;

m = cos(2*pi*fm*t);

figure('Color','w','Position',[100 100 900 600])

subplot(2,1,1)
plot(t,m,'b','LineWidth',2)
hold on
plot(t,-m,'b--')
xline(0,'k','LineWidth',2)
title('Message and Envelope')
ylim([-1.2 1.2]); grid on

subplot(2,1,2)
h_s = animatedline('Color','m','LineWidth',2);
hold on
plot(t,m,'k--')
plot(t,-m,'k--')
xline(0,'k','LineWidth',2)
ylim([-1.2 1.2]); grid on
title('DSB-SC Signal — Observe value at t = 0')

for phi_deg = 0:15:180

    phi = phi_deg*pi/180;
    s = m .* cos(2*pi*fc*t + phi);

    clearpoints(h_s)

    for k = 1:length(t)
        addpoints(h_s,t(k),s(k))
        drawnow
    end

    sgtitle(['Phase \phi = ' num2str(phi_deg) ...
        '°,   s(0) = ' num2str(cos(phi),3)])

    pause(0.7)

end

%% =========================================================
case 8   % Continuous phase change

t=0:0.0005:0.2;
fm=5; fc=50;

figure('Color','w')
h=animatedline('LineWidth',2);
ylim([-1.2 1.2]); grid on
title('Continuous phase change')

for phi=0:0.2:pi
    clearpoints(h)
    m=cos(2*pi*fm*t);
    c=cos(2*pi*fc*t+phi);
    s=m.*c;

    for k=1:length(t)
        addpoints(h,t(k),s(k))
        drawnow
    end
end


%% =========================================================
case 9   % AM vs DSB-SC comparison

dt=0.0005;
t=0:dt:0.2;
fm=5; fc=50;

m=cos(2*pi*fm*t);
c=cos(2*pi*fc*t);

figure('Color','w')

subplot(2,1,1)
h1=animatedline('LineWidth',2); hold on
plot(t,1+m,'k--')
plot(t,-(1+m),'k--')
ylim([-2 2]); grid on
title('AM')

subplot(2,1,2)
h2=animatedline('LineWidth',2); hold on
plot(t,m,'k--')
plot(t,-m,'k--')
ylim([-1.2 1.2]); grid on
title('DSB-SC')

for k=1:length(t)
    addpoints(h1,t(k),(1+m(k))*c(k))
    addpoints(h2,t(k),m(k)*c(k))
    drawnow
end

end