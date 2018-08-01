function mstruct = HBAR_materials(matid, parameters)
 
    match = regexp(matid, '(\w+x)(\d+)?$', 'tokens', 'once');
    if ~isempty(match)
        matid = char(match{1});
        x = char(match{2});
        x = str2double(char(x)) / 100;
    end

    % default values: 
    Poisson      = 0;   % Poisson's Ratio: nu             (nu)
    PiezoCst     = 0;   % Piezoelectric Coefficient: e    (e)
    Resistivity  = inf; % Electrical Resistivity: rho     (r)
    Coupling     = 0;   % Electromechanical Coupling Coefficient: (kt2)
    Permittivity = 1;   % Dielectric Constant: epsilon_r  (eps)
    Losses       = 0;   % Dielectric Losses: tan(delta)   (td)
    
    switch matid
        case {'Si', 'Silicon'}
            Name      = 'Silicon';      %  
            Symbol    = 'Si';           %
            Density   = 2330;           % d
            Qfactor   = 660;            % q
            Stiffness = 166e9;          % c
            Permittivity = 11.68;       % eps
            
        case {'Pt', 'Platinum'}
            Name      = 'Platinum'; 
            Symbol    = 'Pt';
            Density   = 21450;
            Qfactor   = 100;
            Poisson   = 0.38;
            Stiffness = 168e9;
            Permittivity = 1;
            Resistivity  = 1.2e-7;
            
        case {'Mo', 'Molybdenum'}
            Name      = 'Molybdenum';
            Symbol    = 'Mo';
            Density   = 10200;
            Qfactor   = 1000;
            Poisson   = 0.29;
            Stiffness = 270e9;
            Resistivity  = 53.4e-9;
            
        case {'Au', 'Gold'}
            Name        = 'Gold';
            Symbol      = 'Au';
            Density     = 19300;
            Qfactor     = 300;
            Stiffness   = 70e9;
            Poisson     = 0.44;
            Resistivity = 2.2e-8;
            
        case {'Ti' 'Titanium'}
            Name      = 'Titanium';
            Symbol    = 'Ti';
            Density   = 4500;
            Qfactor   = 1000;
            Poisson   = 0.34;
            Stiffness = 115.7e9;
            Resistivity  = 420e-9;
        
        case {'AlN', 'AluminumNitride'}
            Name      = 'Aluminum Nitride';
            Symbol    = 'AlN';
            Density   = 3260;       % rho   [kg/m3]
            Qfactor   = 1000;       % Q     [1]
            Stiffness = 385e9;      % c33   [N/m2]
            PiezoCst  = 1.55;       % e33   [C/m2]
            Permittivity = 10.3;    % eps_r [1]     tanDelta < 0.01 has no effect at all
            Coupling  = 6.5e-2;
            
        case {'ASN15', 'Aluminum-Scandium15-Nitride'}
            Name      = 'Aluminum Scandium Nitride';
            Symbol    = 'ASN';
            Density   = 3280;       % rho   [kg/m3]
            Qfactor   = 1000;       % Q     [1]
            Stiffness = 276.845e9;  % c33   [N/m2]
            PiezoCst  = 1.9;        % e33   [C/m2]
            Permittivity = 10.3;    % eps_r [1]     tanDelta < 0.01 has no effect at all
            Coupling  = 12.3e-2;
			
		case {'ASN26', 'Aluminum-Scandium26-Nitride'}
            Name      = 'Aluminum Scandium Nitride';
            Symbol    = 'ASN';
            Density   = 3295;       % rho   [kg/m3]
            Qfactor   = 1000;       % Q     [1]
            Poisson   = 0;          % nu    [1]
            Stiffness = 300e9;  % c33   [N/m2]
            PiezoCst  = 2.3;        % e33   [C/m2]
            Permittivity = 13.1;    % eps_r [1]     tanDelta < 0.01 has no effect at all
            Resistivity  = inf;     % r     [Ohm*m]
            Coupling  = 16.5e-2;
            
        case {'ASNx'}
            Name      = 'Aluminum Scandium Nitride';
            Symbol    = 'ASN';
            Density   = 3200*(1 + 0.097*x + 0.036*x^2)+30;       % Vladimir's script 
            Qfactor   = 1000;       % Q     [1]
            Poisson   = 0;          % nu    [1]
            %Stiffness = 1e9*(385*(1-x) -23.8*x - 101.4*x*(1-x)); % Caro et al 2015
            %PiezoCst  = 1.46*(1-x) + 8.193*x - 5.912*x*(1-x);    % Caro et al 2015
            Stiffness = 3532*(1-1.305*x+0.207*x^2)*1e8; % Vladimir's script 
            PiezoCst  = 1.483*(1+1.269*x+2.357*x^2);    % Vladimir's script 
            Permittivity = 29.649*x^2 + 15.883*x + 10.515; % our measurements   
            Resistivity  = inf;     % r     [Ohm*m]            
            Losses  = 0.04;
            
            epsf = Permittivity*8.854187817e-12;
            c33E = Stiffness;
            e    = PiezoCst;
            epsS = epsf - e^2/c33E;
            c33D = c33E + e^2/epsS; 
            Coupling = e^2/(c33D*epsS);
            
            
        case {'SiO2', 'SiliconDioxide'}
            Name      = 'Silicon Dioxide';
            Symbol    = 'SiO2';
            Density   = 2197;
            Qfactor   = 500;
            Poisson   = 0.17;
            Stiffness = 74e9;
            PiezoCst  = 0;
            Permittivity = 3.9;
            Resistivity  = inf;
            Coupling  = 0;
    end
    
    mstruct = var2struct(Name, Symbol, Density, Qfactor, Poisson, ...
        Stiffness, PiezoCst, Permittivity, Losses, Resistivity, Coupling );
end
