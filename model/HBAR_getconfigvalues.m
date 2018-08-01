function values = HBAR_getconfigvalues(config, parameters)
    
    layernames = {config.layers.name};
    maternames = {config.materials.Symbol};
    
    proptokens = {'d' 'Q' 'nu' 'c' 'e' 'eps' 'r'};
    propnames = {'Density' 'Qfactor' 'Poisson' 'Stiffness' 'PiezoCst' ... 
        'Permittivity' 'Resistivity'};
    
    tnames = strcat('t', layernames);
    
    nLayers = length(config.layers);
    for k = 1:nLayers
        lcodes{k} = strcat(proptokens, layernames{k});
        mcodes{k} = strcat(proptokens, maternames{k});
    end
    
    if ~iscell(parameters)
        parameters = {parameters};
    end
    
    nParams = length(parameters);
    values  = nan(1, nParams);
    
    for k = 1:nParams; 
        name  = parameters{k};
        
        switch name
            case {'R1' 'R2' 'R3'};
                values(k) = config.geometry.(name);    
            case 'Rc'
                values(k) = config.contact.resistance;
            case 'Lc'
                values(k) = config.contact.inductance;
            case 'xC0'
                values(k) = config.xC0;           
                
            case tnames
                i = find(strcmp(name, tnames), 1);
                values(k) = config.layers(i).thickness;
                
            otherwise
                
                for j = 1:nLayers
                    i = strcmp(lcodes{j}, name) | strcmp(mcodes{j}, name);
                    i = find(i, 1);
                    if ~isempty(i)
                        propname = propnames{i};
                        values(k) = config.materials(j).(propname);
                    end
                end
        end
    end
end