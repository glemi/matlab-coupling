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

Pt.density  = 21450*(1+0.01*1i);
Pt.velocity = 3829;
Pt.young    = 168e9;
Pt.poisson  = 0.38;

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