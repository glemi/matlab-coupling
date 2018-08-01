function config = HBAR_configure(config, parameters)
    
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