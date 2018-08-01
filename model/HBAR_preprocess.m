function hbar = HBAR_preprocess(config)    
    R1 = config.geometry.R1;
    R2 = config.geometry.R2;
    R3 = config.geometry.R3;
    dR = config.geometry.deltaR;
    
    Acntr = pi*(R1+dR)^2;      % Area of central electrode
    Aring = pi*((R3+dR)^2-(R2-dR)^2); % Area of ring electrode
    Aeff  = Acntr*Aring/(Acntr+Aring)*(1+config.xC0);

    TL = struct([]);
    BL = struct([]);
    
    n = length(config.layers);
    for k = 1:n
        
        laydata = config.layers(k);
        matdata = config.materials(k);
        
        material = mkmaterial(matdata);
        layer    = mklayer(laydata, material, Acntr);
        
        switch layer.type(1)
            case 'T',  TL = [TL layer];
            case 'B',  BL = [BL layer];
        end
        switch layer.type
            case 'TE', TE = layer; 
            case 'PL', PL = layer;
            case 'BE', BE = layer;
        end
    end
    
    eps    = PL.material.Permittivity;
    C0     = eps*Aeff/PL.thickness;
    Rc     = config.contact.resistance;
    Lc     = config.contact.inductance;
    Rel    =   TE.material.Resistivity*TE.thickness/Acntr ... % what?
             + BE.material.Resistivity*TE.thickness/Aring; 
    
    hbar.TE    = TE; % Top Electrode
    hbar.PL    = PL; % Piezo Layer
    hbar.BE    = BE; % Bottom Electrode
    hbar.TL    = TL; % Top Layers
    hbar.BL    = BL; % Bottom Layers
    
    hbar.layers= [TL(:); PL(:); BL(:)];
    hbar.Acntr = Acntr;
    hbar.Aring = Aring;
    hbar.Aeff  = Aeff;
    hbar.C0    = C0;
    hbar.Rc    = Rc;
    hbar.Lc    = Lc;
    hbar.Rel   = Rel;
end

function layer = mklayer(laydata, material, area)
    t = laydata.thickness;
    d = material.Density;
    v = material.Velocity;
    
    Zac     = area*d*v;
    phase   = @(f) 2*pi*f*t/v;
    
    function T = tmatrix(f)
        ph = phase(f);
        ph = reshape(squeeze(ph), [1 1 numel(f)]);
        T  = [cos(ph), 1i*Zac*sin(ph); 1i/Zac*sin(ph), cos(ph)];
    end

    layer = laydata;
    layer.area        = area;
    layer.material    = material;
    layer.Zac         = Zac;
    layer.phase       = phase;
    layer.TMatrix     = @tmatrix;
end

function matdata = mkmaterial(matdata)
    epsr = matdata.Permittivity;
    c    = matdata.Stiffness;
    e    = matdata.PiezoCst;
    nu   = matdata.Poisson;
    Q    = matdata.Qfactor;
    rho  = matdata.Density;
    
    % this has no effect if nu = 0;
    c33E = c*(1-nu)/((1+nu)*(1-2*nu)); 
    
    eps0 =  8.854187817e-12;
    epsf = epsr*eps0;
    epsS = epsf - e^2/c33E;
    
    c33D = c33E + e^2/epsS; % applies for AlN only
    c33D = c33D*exp(1i/(2*Q)); % losses 
    
    kt2  = e^2/(c33D*epsS);
    v    = sqrt(c33D/rho);
    
    matdata.Permittivity = epsS;
    matdata.Stiffness = c33D;
    matdata.Velocity = v;
    matdata.Coupling = kt2;
end