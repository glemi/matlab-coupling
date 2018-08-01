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
    config = modify(config, parameters);
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
    R0     = hbar.R0;
    kt2    = real(hbar.PL.material.Coupling);
    
    ZC0    = 1./(2i*pi*f*C0);
    ZL0    = 2i*pi*f*hbar.L0;
    %t = toc; fprintf('electric impedances:  %f\n', t); tic;
    
    N = (ztop+zbot).*sin(gamma) + 2*1i*(1-cos(gamma));
    D = (ztop+zbot).*cos(gamma) + 1i*(1+ztop.*zbot).*sin(gamma);
    Z = ZC0.*(1-kt2./gamma.*N./D) + ZL0 + R0;
    %t = toc; fprintf('total el impedance:   %f\n', t); tic;
    
    data = var2struct(f, gamma, kt2, Ztop, Zbot, C0, ZC0, R0);
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

function config = modify(config, parameters)
    
    layernames = {config.layers.name};
    maternames = {config.materials.Symbol};

    parnames  = parameters(1:2:end-1);
    parvalues = parameters(2:2:end);
    
    proptokens = {'d' 'Q' 'nu' 'c' 'e' 'eps' 'r'};
    propnames = {'Density' 'Qfactor' 'Poisson' 'Stiffness' 'PiezoCst' ... 
        'Permittivity' 'Resistivity'};
    
    tnames = strcat('t', layernames);
    
    nLayers = length(config.layers);
    for k = 1:nLayers
        lcodes{k} = strcat(proptokens, layernames{k});
        mcodes{k} = strcat(proptokens, maternames{k});
    end 
    
    nParams = length(parnames);
    for k = 1:nParams;
        name  = parnames{k};
        value = parvalues{k};
        
        switch name
            case {'R1' 'R2' 'R3'};
                config.geometry.(name) = value;    
            
            case 'R0'
                config.contact.resistance = value;
            
            case 'xC0'
                config.xC0 = value;
                
            case tnames
                i = find(strcmp(name, tnames), 1);
                config.layers(i).thickness = value;
                
            otherwise
                
                for j = 1:nLayers
                    i = strcmp(lcodes{j}, name) | strcmp(mcodes{j}, name);
                    i = find(i, 1);
                    if ~isempty(i)
                        propname = propnames{i};
                        config.materials(j).(propname) = value;
                    end
                end
        end
    end
end