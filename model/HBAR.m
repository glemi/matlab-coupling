%[Zel_HBAR, Zel_FBAR, params] = HBAR(f, varargin)
function [Zel_HBAR, Zel_FBAR, params] = HBAR(f, varargin)
    
    %fprintf('\n\nrstor_model\n');
    
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
    
    assignargs(varargin{:});
    
    %% Frequency 
    f        = reshape(squeeze(f), [1 1 numel(f)]);
    omega    = f*2*pi;
    N        = length(f);
    Zel      = zeros(1,N);
    
    assignargs(varargin{:});
    
    %% Material data
    % quality factor
    
    Qpiezo  = 1000; % AlScN
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
    
    assignargs(varargin{:});
    
    CD33    = (C33+E33^2/Eps33)*exp(1i/(2*Qpiezo)); %Stiffened stiffness 
    
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
    
    %fprintf('kt2 = %.1f%%\n', k2*100);
    assignargs(varargin{:});
    
    %% Velocity in layers
    VPiezo  = (CD33/denPiezo)^0.5;  %Longitudinal velocity in AlScN 
    VSi     = (CSi33/denSi)^0.5;    %Longitudinal velocity in Si 
    VSiO    = ((YoSiO*(1-nuSiO))/((1+nuSiO)*(1-2*nuSiO)*denSiO))^0.5; %Longitudinal velocity in SiO2
    VPt     = ((YoPt*(1-nuPt))/((1+nuPt)*(1-2*nuPt)*denPt))^0.5; %Longitudinal velocity in electrodes (Pt)

    %fprintf('Si resonance:       %.2fMHz\n', real(VSi/tSi/2)/1e6);
    %fprintf('Piezo resonance:    %.2fMHz\n', real(VPiezo/tPiezo/2)/1e6);
    
    
    %% Phase delay beta
    phPiezo  = omega*tPiezo/VPiezo;  %Phase delay AlScN
    phSiO    = omega*tSiO/VSiO;      %Phase delay of SiO2 longitudinal mode
    phSi     = omega*tSi/VSi;        %Phase delay of Si longitudinal mode
    phTopEl  = omega*tTopE/VPt;       %Phase delay of top electrodes longitudinal mode
    phBotEl  = omega*tBotE/VPt;       %Phase delay of bottom electrodes longitudinal mode

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

    %% Matrix and acoustic impedance
    Ztop    = 1i*ZPt*tan(phTopEl); %Top impedance
    ZTOP    = Ztop/ZPiezo; %Top normilized impedance

    TBotEl  = [ cos(phBotEl),            1i*ZPt*sin(phBotEl);...
                1i*ZPt^-1*sin(phBotEl),  cos(phBotEl)         ];
    
    TSiO    = [ cos(phSiO),              1i*ZSiO*sin(phSiO);...
                1i*ZSiO^-1*sin(phSiO),   cos(phSiO)           ];
    
    TSi     = [ cos(phSi),               1i*ZSi*sin(phSi);...
                1i*ZSi^-1*sin(phSi),     cos(phSi)            ];
    
    B    = mmxtimes(mmxtimes(TBotEl,TSiO),TSi);
    Zbot = B(1,2,:)./B(2,2,:);   %Bottom impedance
    ZBOT = Zbot/ZPiezo;     %Bottom normilized impedance
            
    %% Electrical impedance
    AA  = (ZTOP+ZBOT).*sin(phPiezo) + 2*1i*(1-cos(phPiezo));
    BB  = (ZTOP+ZBOT).*cos(phPiezo) + 1i*(1+ZTOP.*ZBOT).*sin(phPiezo);
    
    Zel_HBAR = 1./(1i*omega*C0).*(1-k2./(phPiezo).*AA./BB) + REl;
    Zel_FBAR = 1./(1i*omega*C0).*(1-k2.*tan(phPiezo/2)./(phPiezo/2)) + REl;
    
    f        = squeeze(f);
    Zel      = squeeze(Zel);
    gamma    = squeeze(phPiezo);    
    
    Zel_FBAR = squeeze(Zel_FBAR);
    Zel_HBAR = squeeze(Zel_HBAR);
    
    params.C33         = C33;         % stifness constants AlScN
    params.Eps33       = Eps33;       % Permittivity constants
    params.E33         = E33;         % piezo constants AlScN
    params.CD33        = CD33;        % Stiffened stiffness 
    params.h33         = h33;   
    params.k2          = k2;          % coupling constant
    
    params.denPiezo    = denPiezo;  % AlScN mass density
    params.denSiO      = denSiO;    % Mass density of SiO2
    params.denSi       = denSi;     % Mass density of Si
    params.denPt       = denPt;     % Mass density of electrodes (Pt)
    
    params.tPiezo      = tPiezo;    % AlScN layer thickness
    params.tTopE       = tTopE;     % Top electrode thickness
    params.tBotE       = tBotE;     % Bottom electrode thickness
    params.tSiO        = tSiO;      % SiO2 thickness
    params.tSi         = tSi;       % Si thickness
    
    params.phPiezo    = squeeze(phPiezo);      %Phase delay AlScN
    params.phSiO2     = squeeze(phSiO);        %Phase delay of SiO2 longitudinal mode
    params.phSi       = squeeze(phSi);         %Phase delay of Si longitudinal mode
    params.phTopEl    = squeeze(phTopEl);      %Phase delay of top electrodes longitudinal mode
    params.phBotEl    = squeeze(phBotEl);      %Phase delay of bottom electrodes longitudinal mode

    params.ZPiezo     = ZPiezo;       %Longitudinal wave impedance of active layer 
    params.ZSiO       = ZSiO;         %Longitudinal wave impedance of SiO2 
    params.ZSi        = ZSi;          %Longitudinal wave impedance of Si
    params.ZPt        = ZPt;          %Longitudinal wave impedance of Pt
    
    params.Ztop       = squeeze(Ztop);
    params.Zbot       = squeeze(Zbot);
    
    params.numerator = squeeze(AA);
    params.denominator = squeeze(BB);

    params.VPiezo  	  = VPiezo;       %Longitudinal velocity in AlScN 
    params.VSi     	  = VSi;          %Longitudinal velocity in Si 
    params.VSiO    	  = VSiO;         %Longitudinal velocity in SiO2
    params.VPt     	  = VPt;          %Longitudinal velocity in Pt
    
    params.frPiezo    = real(VPiezo/tPiezo/2);
    params.frSi       = real(VSi/tSi/2);
    params.frSiO2     = real(VSiO/tSiO/2);
    params.frTopEl    = real(VPt/tTopE/2);
    params.frBotEl    = real(VPt/tBotE/2);
    
    params.gamma      = gamma;

    
%     fig rstordebug:HBAR; plot(f, abs(Zel_HBAR)); xscale log; yscale log;
%     fig rstordebug:FBAR; plot(f, abs(Zel_FBAR)); xscale log; yscale log;
end
