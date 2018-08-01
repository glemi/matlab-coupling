% function [Z, data] = HBAR_v3(f, config, varargin)
%  
% HBAR Model v3
%
function [Z, data] = HBAR_v3(f, config, varargin)
    opt config char|struct HBAR_config.txt;

    parameters = varargin;
    f = [f(:)]; %#ok
    
    config = getconfig(config);
    config = HBAR_parameters(config, parameters);
    %t = toc; fprintf('load & modify config: %f\n', t); tic;
    
    hbar   = HBAR_preprocess(config);
    %t = toc; fprintf('preprocess:           %f\n', t); tic;
    
    ZPiezo = hbar.PL.Zac; 
    gamma  = hbar.PL.phase(f);
    Ztop   = zcompute(f, hbar.TL);
    Zbot   = zcompute(f, hbar.BL);
    ztop   = Ztop/ZPiezo;
    zbot   = Zbot/ZPiezo;
    %t = toc; fprintf('compute impedances:   %f\n', t); tic;
    
    C0     = hbar.C0;
    Rc     = hbar.Rc;
    Lc     = hbar.Lc;
    kt2    = (hbar.PL.material.Coupling);
    
    ZC0    = 1./(2i*pi*f*C0);
    %ZLc    = 2i*pi*f*Lc;
    Z0     = ZC0; %.*ZLc./(ZC0+ZLc) +Rc;
    %t = toc; fprintf('electric impedances:  %f\n', t); tic;
    
    N = (ztop+zbot).*sin(gamma) + 2*1i*(1-cos(gamma));
    D = (ztop+zbot).*cos(gamma) + 1i*(1+ztop.*zbot).*sin(gamma);
    H = (1 - kt2./gamma.*N./D);
    Z = Z0.*H;
    
    %t = toc; fprintf('total el impedance:   %f\n', t); tic;
    
    data = var2struct(f, gamma, kt2, Z, H, Ztop, Zbot, C0, ZC0, Rc);
    data.config = config;
    data.layers = hbar.layers;
    %t = toc; fprintf('prepare output data:  %f\n', t); tic;
end

function z = zcompute(f, layers)
    n = length(layers);
    
    T = layers(1).TMatrix(f);
    for k = 2:n
        T1 = layers(k).TMatrix(f);
        T  = mmxtimes(T, T1);
    end
    
    z = squeeze(T(1,2,:)./T(2,2,:));
end

function config = getconfig(config)
    persistent lastfile lastconfig;
    
    if ischar(config) && exist(config, 'file')
        filename = config;
        if strcmp(filename, lastfile)
            config = lastconfig;
        else
            config = HBAR_loadconfig(filename);
            lastfile = filename;
            lastconfig = config;
        end
    elseif ~isstruct(config)
        error 'invalid config data or file path';
    end
end