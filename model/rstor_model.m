function [Zel, f, gamma, k2] = rstor_model(varargin)
    
    fprintf('\n\nrstor_model\n');
    
    assignargs(varargin{:});

    %% Geometry
    R1 = 95e-6;         %Electrode radius
    R2 = 210e-6;
    R3 = 300e-6;
    S  = pi*R1^2;      % Area of central electrode
    S2 = pi*(R3-R2)^2; % Area of ring electrode
    
    %% Thicknesses
    tPiezo  = 800e-9;  % AlScN layer thickness
    tTopE   = 100e-9;   % Top electrode thickness
    tBotE   = 200e-9;  % Bottom electrode thickness
    tSiO    = 200e-9;  % SiO2 thickness
    tSi     = 725e-6;  % Si thickness
    
    Eps0    = 8.854187817e-12; %vacuum permittivity

    RhoPt   = 1.2e-7;  % Pt resistance Ohm*m
    REl     = RhoPt*tTopE/S+RhoPt*tBotE/S2; %Pt electrodes resistivity
    
    fstart = 1e9;
    fstep = 100e3;
    fstop = 10e9;
    
    assignargs(varargin{:});
    
    %% Frequency 
    f(1,1,:) = fstart:fstep:fstop;
    omega    = f*2*pi;
    N        = length(f);
    Zel      = zeros(1,N);
    
    assignargs(varargin{:});
    
    %% Material data
    % quality factor
    
    QPiezo  = 1000; % AlScN
    QSiO    = 500;  % SiO2
    QSi     = 660;  % Si
    QPt     = 100;  % Pt
    
    % QPiezo = 1500;
    % QSiO2  = 500;
    % QSi    = 3000;
    % QPt    = 200;
    
    tanDelta = 0.002;
    
    
    assignargs(varargin{:});
    
    denPiezo = 3280; % AlScN mass density
    denSiO  = 2197;  % Mass density of SiO2
    denSi   = 2330;  % Mass density of Si
    denPt   = 21450; % Mass density of electrodes (Pt)

    YoSiO   = 74e9*exp(1i/(2*QSiO)); %Young's modulus of SiO2
    YoPt    = 168e9*exp(1i/(2*QPt)); %Young's modulus of electrodes (Pt) 
    %YoSiO   = 70e9; %Young's modulus of SiO2
    %YoPt    = 168e9; %Young's modulus of electrodes (Pt) 

    nuSiO   = 0.17; %Poisson's ratio of SiO2
    nuPt    = 0.38; %Poisson's ratio of electrodes (Pt) 

    CSi33   = 166e9*(1+tan(1i/(2*QSi))); % stifness constants Si
    %CSi33   = 166e9; % stifness constants Si

    C33     = 276.845e9;                    % stifness constants AlScN
    Eps33   = 14*Eps0*exp(-1i*tanDelta);    % Permittivity constants
    E33     = 1.9;                          % 1.883; % piezo constants AlScN
    CD33    = (C33+E33^2/Eps33)*exp(1i/(2*QPiezo)); %Stiffened stiffness 
    
    %C33     = 276.845e9; % stifness constants AlScN
    %E33     = 1.883; % piezo constants AlScN
    %Eps33   = 14.1*Eps0; %Permittivity constants
    %CD33    = C33+E33^2/Eps33; %Stiffened stiffness 

    assignargs(varargin{:});
    
    %% Coupling square k2 and another parameters
    h33      = E33/Eps33; 
    k2       = real(E33^2/(CD33*Eps33));
    %h33     = E33/Eps33; 
    %k2      = h33^2*Eps33/CD33;
    
    fprintf('kt2 = %.1f%%\n', k2*100);
    assignargs(varargin{:});
    
    %% Velocity in layers
    VPiezo  = (CD33/denPiezo)^0.5;  %Longitudinal velocity in AlScN 
    VSi     = (CSi33/denSi)^0.5;    %Longitudinal velocity in Si 
    VSiO    = ((YoSiO*(1-nuSiO))/((1+nuSiO)*(1-2*nuSiO)*denSiO))^0.5; %Longitudinal velocity in SiO2
    VPt     = ((YoPt*(1-nuPt))/((1+nuPt)*(1-2*nuPt)*denPt))^0.5; %Longitudinal velocity in electrodes (Pt)

    fprintf('Si resonance:       %.2fMHz\n', real(VSi/tSi)/1e6);
    fprintf('Piezo resonance:    %.2fMHz\n', real(VPiezo/tPiezo)/1e6);
    
    
    %% Phase delay beta
    
    betaPiezo  = omega*tPiezo/VPiezo;  %Phase delay AlScN
    betaSiO    = omega*tSiO/VSiO;      %Phase delay of SiO2 longitudinal mode
    betaSi     = omega*tSi/VSi;        %Phase delay of Si longitudinal mode
    betaTopEl  = omega*tTopE/VPt;       %Phase delay of top electrodes longitudinal mode
    betaBotEl  = omega*tBotE/VPt;       %Phase delay of bottom electrodes longitudinal mode

    fcPiezo = 1/2*VPiezo/tPiezo;
    fcOxide = 1/2*VSiO/tSiO;
    fcSubst = 1/2*VSi/tSi;
    fcTopEl = 1/2*VPt/tTopE;
    fcBotEl = 1/2*VPt/tBotE;
    
    
    %% Acoustic impedance of layers Z
    ZPiezo  = S*denPiezo*VPiezo; %Longitudinal wave impedance of active layer 
    ZSiO    = S*denSiO*VSiO;       %Longitudinal wave impedance of SiO2 
    ZSi     = S*denSi*VSi;          %Longitudinal wave impedance of Si
    ZPt     = S*denPt*VPt;          %Longitudinal wave impedance of Pt
    
    assignargs(varargin{:});
    
    C1 = Eps33*S/tPiezo;  %Capacitance of the central electrode
    C2 = Eps33*S2/tPiezo; %Capacitance of the side electrode
    C0 = C1*C2/(C1+C2);   %Series connection
    %C0 = 3e-12;
    %Sigma = h33^2*(omega*ZPiezo).^-1*C0;
    
    %sigma   = k2./betaPiezo;
    %C0      = Eps33*S/tPiezo; %AlScN layer capacitance

    %% Matrix and acoustic impedance
    Ztop    = 1i*ZPt*tan(betaTopEl); %Top impedance
    ZTOP    = Ztop/ZPiezo; %Top normilized impedance

    TBotEl  = [ cos(betaBotEl),            1i*ZPt*sin(betaBotEl);...
                1i*ZPt^-1*sin(betaBotEl),  cos(betaBotEl)         ];
    
    TSiO    = [ cos(betaSiO),              1i*ZSiO*sin(betaSiO);...
                1i*ZSiO^-1*sin(betaSiO),   cos(betaSiO)           ];
    
    TSi     = [ cos(betaSi),               1i*ZSi*sin(betaSi);...
                1i*ZSi^-1*sin(betaSi),     cos(betaSi)            ];
    
    B    = mmxtimes(mmxtimes(TBotEl,TSiO),TSi);
    Zbot = B(1,2,:)./B(2,2,:);   %Bottom impedance
    ZBOT = Zbot/ZPiezo;     %Bottom normilized impedance
            
%     A = [ cos(betaPiezo), -1i*ZPiezo*sin(betaPiezo), -1i*(h33./omega).*(1-cos(betaPiezo));...   %Piezoelectric matrix
%          -1i*ZPiezo^-1*sin(betaPiezo), cos(betaPiezo),(h33./omega)*ZPiezo^-1.*sin(betaPiezo);...
%          (h33./omega)*ZPiezo^-1.*sin(betaPiezo), -1i*(h33./omega).*(1-cos(betaPiezo)), -1i*(1-sigma.*sin(betaPiezo))./(omega*C0)];
% 
%     Ztop =  1i*ZPt*tan(betaTopEl); %Top impedance
% 
%     Zbot = 1i*ZPt*...
%          (ZSi*(ZPt-ZSiO*cot(betaBotEl).*cot(betaSiO)) - ZSiO*cot(betaSi).*(ZSiO*cot(betaBotEl)+ZPt*cot(betaSiO)))...
%        ./(ZSi*(ZPt*cot(betaBotEl)+ZSiO*cot(betaSiO))  + ZSiO*cot(betaSi).*(ZSiO-ZPt*cot(betaBotEl).*cot(betaSiO))); %Bottom impedance

    %% Electrical impedance
    %Zel = (A(3,1,:).*(A(1,3,:)-A(2,3,:).*Zbot).*(Ztop+A(3,2,:)))./(Zbot.*A(2,1,:).*(Ztop+A(2,2,:))-A(1,1,:).*(Ztop+A(1,2,:)))+A(3,3,:);
    AA  = (ZTOP+ZBOT).*sin(betaPiezo) + 2*1i*(1-cos(betaPiezo));
    BB  = (ZTOP+ZBOT).*cos(betaPiezo) + 1i*(1+ZTOP.*ZBOT).*sin(betaPiezo);
    Zel = 1./(1i*omega*C0).*(1-k2./(betaPiezo).*AA./BB) + REl;
    
    Zel_FBAR = 1./(1i*omega*C0).*(1-k2.*tan(betaPiezo/2)./(betaPiezo/2)) + REl;
    Zel_HBAR = Zel;
    
    f = squeeze(f);
    Zel = squeeze(Zel);
    gamma = squeeze(betaPiezo);    
    
    Zel_FBAR = squeeze(Zel_FBAR);
    Zel_HBAR = squeeze(Zel_HBAR);
%     
%     fig rstordebug:HBAR; plot(f, abs(Zel_HBAR)); xscale log; yscale log;
%     fig rstordebug:FBAR; plot(f, abs(Zel_FBAR)); xscale log; yscale log;
%     
%     %pattern = 'R = %.2e Ohm \t Q = %.0f  \t C0 = %.0f pF  \t fr = %.1fGHz  \t k = %.3f\n';
%     %fprintf(pattern, R, Q, C0*1e12, fr*1e-9, k);
%     fprintf('Resonator Model:\n');
%     fprintf(' C0 = %.0f pF\n', C0*1e12);
%     fprintf(' k2 = %.3f\n', k2);    

    fig rstordebug:HBAR; clf;
    subplot(2,3,1); title Z_{bottom}
    plot(f, abs(squeeze(ZBOT))); xscale log; yscale log;
    subplot(2,3,2); title Z_{top}
    plot(f, abs(squeeze(ZTOP))); xscale log; yscale log;
    subplot(2,3,3); 
    plot(f, abs(squeeze((ZTOP+ZBOT)./(1 + ZTOP.*ZBOT)))); xscale log; yscale log;
end

