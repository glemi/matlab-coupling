clear;
set(0, 'defaultaxesnextplot', 'add');
set(0, 'defaultaxesbox', 'on');
set(0, 'defaultaxesxscale', 'log');
set(0, 'defaultaxesyscale', 'log');

addpath ../universal;

%% frequency range
n = 10000;
f = logspace(6, 11, n);
w = 2*pi*f;

%% Other parameters
C0 = 250e-12; % Static capacitance: 250pF
kt = 1;       % Coupling factor: ?

%% material properties;
Si.density = 2328.3;
Si.velocity = 8433.2;

Pt.density = 21450;
Pt.velocity = 3829;

AlScN.density = 3280;
AlScN.velocity = 9721;

%% restonator stack

r = 500e-6;
A = r^2*pi;

% top electrode
TE.thickness    = 100e-9;
TE.Z            = A*Pt.density*Pt.velocity;
TE.phaseDelay   = w*TE.thickness/Pt.velocity;

% piezo layer
PL.thickness    = 1000e-9;
PL.Z            = A*AlScN.density*AlScN.velocity;
PL.phaseDelay   = w*PL.thickness/AlScN.velocity;

% bottom electrode
BE.thickness    = 100e-9;
BE.Z            = A*Pt.density*Pt.velocity;
BE.phaseDelay   = w*BE.thickness/Pt.velocity;

% substrate
Sb.thickness    = 525e-6;
Sb.Z            = A*Si.density*Si.velocity;
Sb.phaseDelay   = w*Sb.thickness/Si.velocity;



%% Compute Input Impedance

Z1 = 1i*TE.Z*tan(TE.phaseDelay);

Z2 = 1i*(Sb.Z*tan(Sb.phaseDelay) + TE.Z*tan(TE.phaseDelay)) ...
      ./ (1 - (Sb.Z/BE.Z)*tan(BE.phaseDelay).*tan(Sb.phaseDelay));

ZC0 = 1./(1i*w*C0);

z1 = Z1/Sb.Z;
z2 = Z2/Sb.Z;


R =  ( (z1 + z2).*sin(PL.phaseDelay) + 2i*(1-cos(PL.phaseDelay))) ...
   ./ (z1 + z2).*cos(PL.phaseDelay) + 1i*(1+z1.*z2).*sin(PL.phaseDelay);

R1 =   sin(PL.phaseDelay) + 2i*(1-cos(PL.phaseDelay)) ...
   ./ cos(PL.phaseDelay) + 1i*sin(PL.phaseDelay);

Zin = ZC0.*(1 - kt^2./PL.phaseDelay.*R);


%% plot

fig rsm:gak; clf; 
plot(f, abs(Zin));
plot(f, abs(ZC0));
%xscale log; yscale log;
axmenu;

fig rsm:gak1; clf; 
plot(f, abs(z1));
title |z_1|
%plot(f, abs(z2));
%xscale log; yscale log;
axmenu;


fig rsm:gak2; clf; 
plot(f, abs(z2));
title |z_2|
%plot(f, abs(z2));
%xscale log; yscale log;
axmenu;

fig rsm:gak3; clf;
plot(f, abs(R1));
%yscale log;


fig rsm:gak4; clf;
findpeaks(abs(Zin-ZC0), f);
[~, fr] = findpeaks(abs(Zin-ZC0), f);


fig rsm:gak5; clf;
plot(fr(2:end), diff(fr), '.');
xUnitTicks 'Hz'

