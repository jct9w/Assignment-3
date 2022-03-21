clear all
close all

% constants %
T = 300;
K = 1.3806504e-23;
Mo = 9.10938215e-31;
Mn = 0.26*Mo;
Vth = sqrt((K*T)/(Mn));
Tmn = 0.2e-12;
%MFP = Tmn * Vth;
q = 1.60217653e-19;


% size of modeling %
L = 200e-9;
W = 100e-9;

% basic modifiable parameters %
esize = 3000; % electron numbers
tstep = 15e-15; % time step
simlen = 50; % simulation timesteps

% initiate electron specific data %
Px = rand(1,esize) * 2e-7;
Py = rand(1,esize) * 1e-7;
edir = rand(1,esize) * 2 * pi;
Vx = Vth * cos(edir);
Vy = Vth * sin(edir);
edata = zeros(4,esize);
edata(1,:) = Px;
edata(2,:) = Py;
edata(3,:) = Vx;
edata(4,:) = Vy;

scatdist = makedist('Normal', 'mu', 0, 'sigma', sqrt(K*T/Mn));

% misc (plot color data) %
carr = rand(1,esize);

% MFP specific variables %
MFP = zeros(1,esize);
MFPx = Px; %temp var
MFPy = Py; %temp var
MTBC = zeros(1,esize);
Pscat = 1 - exp(-tstep/Tmn);
Tsum = 0; %temp var

%% Assignment 3 p1

Voltx = 0.1;
Volty = 0;

eEx = Voltx/L;
eEy = Volty/W;
eFx = eEx*q;
eFy = eEy*q;

Acc = zeros(2,esize);
Acc(1,:) = eFx/Mn;
Acc(2,:) = eFy/Mn;

ec = 10e15;
edcd = zeros(1,simlen);
tdenvar = zeros(1,esize);

%%

% main function %
for n = 1 : simlen

    temp = sqrt(edata(3,:).^2 + edata(4,:).^2);

    % added dx and dy from EMModeling for faster access %
    dx = edata(3,:) * tstep;
    dy = edata(4,:) * tstep;

    % Assignment 3 acceleration
    if Acc(1,:) ~= 0
        edata(3,:) = edata(3,:) + Acc(1,:)*tstep;
    end
    if Acc(2,:) ~= 0
        edata(4,:) = edata(4,:) + Acc(1,:)*tstep;
    end
    
    % probablity of scattering %
    scat = rand(1,esize) < Pscat;
    edata(3:4,scat) = random(scatdist,[sum(scat),2])';

    % scattered temp and avg temp calculation %
    nVth = sqrt(edata(3,:).^2 + edata(4,:).^2);
    avgV = sum(nVth) / esize;
    nT = T + ((Mn * (avgV.^2)) / K / esize / 2);
    Tsum = Tsum + nT;
    avgT = Tsum / n;

    % MFP and MTBC %
    MFP(scat) = MFP(scat) + sqrt((dx(scat)-MFPx(scat)).^2+(dy(scat)-MFPy(scat)).^2);
    MTBC(scat) = MTBC(scat) + sqrt ((dx(scat)-MFPx(scat)).^2+(dy(scat)-MFPy(scat)).^2)./nVth(scat);

    % new cordinates calculation %
    nPx = edata(1,:) + dx;
    nPy = edata(2,:) + dy;
    
    % check if new cordinates excedes boundary %
    uplim = nPy > W;
    rtlim = nPx > L;
    lflim = nPy < 0;
    dwlim = nPx < 0;
    % actions taken %
    nPy(uplim) = 2 * W - nPy(uplim);
    edata(4,uplim) = -1 * edata(4,uplim);
    nPx(rtlim) = nPx(rtlim) - L;
    nPy(lflim) = -1 * nPy(lflim);
    edata(4,lflim) = -1 * edata(4,lflim);
    nPx(dwlim) = nPx(dwlim) + L;

    % 2d traj plot %
    figure(1)
    scatter(nPx,nPy,2,carr)
    title '2-D plot of particle trajectories'
    xlabel 'x'
    ylabel 'y'
    xlim([0 L]);
    ylim([0 W]);
    hold on
    grid on

    % temp plot %
%     figure(2)
%     subplot(2,1,1)
%     scatter(n,nT,2)
%     title 'Temperature plot'
%     hold on
%     xlim([0 simlen])
%     ylim([300 330])
% 
%     subplot(2,1,2)
%     scatter(n,avgT,2)
%     title 'average temperature trend plot'
%     hold on
%     xlim([0 simlen])
%     ylim([300 330])

    % cordinates update %
    edata(1,:) = nPx;
    edata(2,:) = nPy;

    % Assignment 3 current density
    Vdrift = mean(temp)*eEx;
    mu = Vdrift/esize;
    edcd(n) =q*ec*mu*W*L

    tdenvar = Mn*(temp.^2)/(K*2);
end

%avgT
%avgMFP = mean(MFP)
%avgMTBC = mean (MTBC)

figure(3)
plot(linspace(0,simlen,simlen),edcd)
title('Current plot')
xlabel('time step')
ylabel('Current')
hold on

% Density
edenmap = [edata(1,:)', edata(2,:)'];
figure(4)
surf(hist3(edenmap,[20,10]))
title 'Electron Density Map';
zlabel 'Number of Electrons per Grid Point';
ylabel 'Y Coordinate (*10^-7m)';
xlabel 'X coordinate (*10^-7m)';

figure(5)
xplot = linspace(0,200e-9,40);
yplot = linspace(0,100e-9,40);
[X, Y] = ndgrid(xplot', yplot');
F = scatteredInterpolant(edata(1,:)', edata(2,:)', tdenvar');
Z = F(X,Y);
surf(X,Y,Z);
title 'Temperature Map';
zlabel 'Temperature';
ylabel 'Y Coordinate (*10^-7m)';
xlabel 'X coordinate (*10^-7m)';


%     figure(4)
%     subplot(2,1,1)
%     histogram(edata(3,:),simlen/25)
%     title 'histogram of vx'
%     subplot(2,1,2)
%     histogram(edata(4,:),simlen/25)
%     title 'histogram of vy'
