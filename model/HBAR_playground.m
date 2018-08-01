function [Z, data] = HBAR_playground(f, varargin)
    parameters = varargin;
    f = [f(:)]; %#ok
    
    R0 = 95e-6;    R1 = 210e-6;    R2 = 300e-6;
    Acntr = pi*R0^2;      % Area of central electrode
    Aring = pi*(R2-R1)^2; % Area of ring electrode

    layer('TopEl', 'Pt',   100e-9, Acntr, parameters);
    layer('Piezo', 'AlN',  800e-9, Acntr, parameters);
    layer('BotEl', 'Pt',   200e-9, Acntr, parameters);
    layer('Oxide', 'SiO2', 200e-9, Acntr, parameters);
    layer('Subst', 'Si',   725e-6, Acntr, parameters);
    
    toplayers   = [TopEl];
    botmlayers  = [BotEl Oxide Subst];
    
    ZPiezo = Piezo.Zac;
    gamma  = Piezo.phase(f);
    Ztop   = zcompute(f, toplayers);
    Zbot   = zcompute(f, botmlayers);
    ztop   = Ztop/ZPiezo;
    zbot   = Zbot/ZPiezo;
    
    eps    = Piezo.material.Permittivity;
    Ccntr  = eps*Acntr/Piezo.thickness; % central electrode
    Cring  = eps*Aring/Piezo.thickness; % ring electrode
    C0     = Ccntr*Cring/(Ccntr+Cring);   %Series connection

    ZC0    = 1./(2i*pi*f*C0);
    Rel    = TopEl.material.Resistivity*TopEl.thickness/Acntr ... % what?
            +BotEl.material.Resistivity*BotEl.thickness/Aring;    % is this right?
    
    kt2    = real(Piezo.material.Coupling);
        
    N = (ztop+zbot).*sin(gamma) + 2*1i*(1-cos(gamma));
    D = (ztop+zbot).*cos(gamma) + 1i*(1+ztop.*zbot).*sin(gamma);
    Z = ZC0.*(1-kt2./gamma.*N./D) + Rel;
    
    data = var2struct(f, gamma, kt2, Ztop, Zbot, Cring, ZC0, Rel);
    data.layers = [TopEl Piezo BotEl Oxide Subst];
end

function z = zcompute(f, layers)
    n = length(layers);
    T = repmat(eye(2), 1, 1, length(f));
 
    for k = 1:n
        T1 = layers(k).TMatrix(f);
        T  = mmxtimes(T, T1);
    end
    
    z = squeeze(T(1,2,:)./T(2,2,:));
end

function lstruct = layer(layid, matid, thickness, area, parameters)
    material = materials(matid, parameters);
    override(layid, parameters);
    
    d = material.Density;
    v = material.Velocity;
    
    Zac     = area*d*v;
    phase   = @(f) 2*pi*f*thickness/v;
    TMatrix = @tmatrix;
    
    function T = tmatrix(f)
        ph = phase(f);
        ph = reshape(squeeze(ph), [1 1 numel(f)]);
        T  = [cos(ph), 1i*Zac*sin(ph); 1i/Zac*sin(ph), cos(ph)];
    end

    lstruct = var2struct(Zac, phase, TMatrix, thickness, area, material);
    if nargout == 0
        assignin('caller', layid, lstruct);
    end
end

function mstruct = materials(matid, parameters)
 
    switch matid
        case {'Si', 'Silicon'}
            Name      = 'Silicon';      %  
            Symbol    = 'Si';           %
            Density   = 2330;           % d
            Qfactor   = 660;            % q
            Poisson   = 0;              % nu
            Stiffness = 166e9;          % c
            PiezoCst  = 0;              % e
            Permittivity = 11.68;       % eps
            Resistivity  = inf;         % r
            
        case {'Pt', 'Platinum'}
            Name      = 'Platinum';
            Symbol    = 'Pt';
            Density   = 21450;
            Qfactor   = 100;
            Poisson   = 0.38;
            Stiffness = 168e9;
            PiezoCst  = 0;
            Permittivity = 1;
            Resistivity  = 1.2e-7;
            
        case {'Mo', 'Molybdenum'}
            Name      = 'Molybdenum';
            Symbol    = 'Mo';
            Density   = 10200;
            Qfactor   = 1000;
            Poisson   = 0.29;
            Stiffness = 270e9;
            PiezoCst  = 0;
            Permittivity = 1;
            Resistivity  = 53.4e-9;
            
        case {'AlN', 'AluminumNitride'}
            Name      = 'Aluminum Nitride';
            Symbol    = 'AlN';
            Density   = 3280;       % rho   [kg/m3]
            Qfactor   = 1000;       % Q     [1]
            Poisson   = 0;          % nu    [1]
            Stiffness = 276.845e9;  % c33   [N/m2]
            PiezoCst  = 1.9;        % e33   [C/m2]
            Permittivity = 14;      % eps_r [1]     tanDelta < 0.01 has no effect at all
            Resistivity  = inf;     % r     [Ohm*m]
            
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
    end
    
    eps0 =  8.854187817e-12;
    
    override(matid, parameters);
    Permittivity = Permittivity*eps0;
    Stiffness = Stiffness*(1-Poisson)/((1+Poisson)*(1-2*Poisson)); % no effect if Poisson = 0;
    Stiffness = Stiffness + PiezoCst^2/Permittivity; % applies for AlN only
    Stiffness = Stiffness*exp(1i/(2*Qfactor)); % losses 
    Velocity  = sqrt(Stiffness/Density);
    
    Coupling  = PiezoCst^2/(Stiffness*Permittivity);
    
    mstruct = var2struct(Name, Symbol, Density, Qfactor, Poisson, ...
        Stiffness, PiezoCst, Permittivity, Resistivity, Velocity, ...
        Coupling);
end

function count = override(id, parameters)
    count = 0;
    short = {'d' 'Q' 'nu' 'c' 'e' 'eps' 'r', 't'};
    vname = {'Density' 'Qfactor' 'Poisson' 'Stiffness' 'PiezoCst' ... 
        'Permittivity' 'Resistivity', 'thickness'};

    n = length(parameters);
    for k = 1:2:n-1;
        aname = parameters{k};
        value = parameters{k+1};
        
        names = strcat(short, id);
        i = find(strcmp(aname, names));
        
        if ~isempty(i)
            assignin('caller', vname{i}, value);
            count = count + 1;
        end
    end
    
end
