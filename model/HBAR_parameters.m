% function output = HBAR_parameters(config, params)
% 
% Read or modify parameters of a HBAR configuration structure. 
%
% usage:
%           % write:
%           params = {'tPiezo' 800e-9 'ePiezo' 1.8};
%           config = HBAR_parameters(config, params);
%           
%           % read:
%           param_names = {'tPiezo' 'cPiezo' 'tTopEl'};
%           values = HBAR_parameters(config, param_names);
%           
% valid parameter names:
%
% Rc            : Contact Resistance
% Lc            : Contact Inductance
%
% R1, R2, R3    : Electrode Radii
% 
% eps<LayerName> or eps<MaterialName> : Permittivity (relative) 
%  td<LayerName> or  td<MaterialName> : Dielectric Losses (tan delta) [1]
%  nu<LayerName> or  nu<MaterialName> : Poisson's Ratio      [1]
%   d<LayerName> or   d<MaterialName> : Density              [kg/m3]
%   Q<LayerName> or   Q<MaterialName> : Quality Factor       [1]
%   c<LayerName> or   c<MaterialName> : Stiffness (c33)      [Pa]
%   e<LayerName> or   e<MaterialName> : Peizo Constant (e33) [1]
%   r<LayerName> or   r<MaterialName> : Resistivity          [Ohm*m]
%   t<LayerName> or   t<MaterialName> : Thickness            [m]
%
% see also: HBAR_loadconfig, HBAR_v3, HBAR_preprocess, HBAR_materials;
function output = HBAR_parameters(config, params)
    if ischar(params)
        params = cellstr(params);
    end
    
    if isempty(params)
        output = config;
        return;
    elseif iscellstr(params)
        parnames = params;
        parvalues = [];
        writemode = false;
    else
        parnames  = params(1:2:end-1);
        parvalues = params(2:2:end);
        writemode = true;
    end

    
    n = length(parnames);
    if writemode
        for k = 1:n
            config = writevalue(config, parnames{k}, parvalues{k});
        end
        output = config;
    else
        m = length(config);
        for j = 1:m
            for k = 1:n
                try 
                    parvalues(j,k) = readvalue(config(j), parnames{k});
                catch 
                    parvalues(j,k) = NaN;
                end
            end
        end
        output = parvalues;
    end
end

function value = readvalue(config, parcode)
    [parname, groupname, layname] = dissect(parcode);
    i = find(getlaynumber(config, layname), 1);
    
    switch groupname
        case 'layers'
            value = config.layers(i).(parname);
        case 'materials'
            value = config.materials(i).(parname);
        case ''
            value = config.(parname);
        otherwise
            value = config.(groupname).(parname);
    end
end

function config = writevalue(config, parcode, value)
    [parname, groupname, layname] = dissect(parcode);
    i = getlaynumber(config, layname);

    switch groupname
        case 'layers'
            [config.layers(i).(parname)] = deal(value);
        case 'materials'
            [config.materials(i).(parname)] = deal(value);
        case ''
            config.(parname) = value;
        otherwise
            config.(groupname).(parname) = value;
    end
    
    switch parcode
        case 'kt2', config.override.kt2 = value;
        case 'C0',  config.override.C0 = value;
    end
end

function ilayer = getlaynumber(config, layname)
    layers = config.layers;
    lnames = {layers.name};
    ltypes = {layers.type};
    mnames = {layers.material};
    
    ilayer = strcmp(layname, lnames);
    iltype = strcmp(layname, ltypes);
    imater = strcmp(layname, mnames);
    ilayer = ilayer | iltype | imater;
    
    if ~any(ilayer) && ~isempty(layname)
        error 'Invalid Layer or Material Code';
    end
end

function [parname, groupname, layname] = dissect(code)
    %% Special case: coupling coefficient
    if strcmp(code, 'kt2')
        parname = 'Coupling';
        groupname = 'materials';
        layname = 'PL';
        return;
    end

    %% general parameters
    parcodes = {'Rc' 'Lc' 'R1' 'R2' 'R3' 'xC0' 'C0'} ;
    parnames = {'resistance', 'inductance', 'R1' 'R2' 'R3' 'xC0' 'C0'};
    groupnames = {'contact' 'contact' 'geometry' 'geometry' 'geometry' '' ''};
    i = strcmp(code, parcodes);
    groupname = [groupnames{i}];
    parname = [parnames{i}];
    layname = '';
    
    if any(i)
        return; 
    end

    %% per-layer parameters
    parcodes = {'eps' 'td' 'nu' 'd' 'Q' 'c' 'e' 'r' 't'};
    parnames = {'Permittivity' 'Losses' 'Poisson' 'Density' 'Qfactor' ...
        'Stiffness' 'PiezoCst' 'Resistivity' 'thickness'};
    
    
    n = length(parcodes);
    for k = 1:n
        parcode = parcodes{k};
        m = length(parcode);
        if strncmp(code, parcode, m)
            layname = code(m+1:end);
            parname = parnames{k};
            break;
        end
    end
    
    switch parcode
        case 't',  groupname = 'layers';
        otherwise, groupname = 'materials';
    end
    
    if isempty(layname)
        error 'invalid parameter code'; 
    end
end