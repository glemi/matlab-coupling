function config = HBAR_loadconfig(filepath, varargin)    
    if nargin < 1 || isempty(filepath)
        filepath = 'HBAR_config.txt';
    end
    
    [config, parameters] = readfile(filepath);
    
    if config.contact.inductance == 0
        config.contact.inductance = 1e50;
    end
    
    [~, filetitle] = fileparts(filepath);
    config.file = filetitle;
    config.filepath = filepath;
    
    config.override.C0 = [];
    config.override.kt2 = [];
    
    config.C0 = 1.7e-12;
    config.xC0 = 0.15;
    
    ubound = config;
    lbound = config;
    
    n = length(parameters);
    for k = 1:n
        par = parameters(k);
        config = HBAR_parameters(config, {par.name par.value});
        ubound = HBAR_parameters(ubound, {par.name par.ubound});
        lbound = HBAR_parameters(lbound, {par.name par.lbound});
    end
    config.ubound = ubound;
    config.lbound = lbound;
    
    config = HBAR_parameters(config, varargin);
end

function [config, parameters] = readfile(filename)
    k = 0; j = 0;
    
    config = struct;
    %parameters = struct([]);
    
    fileID = fopen(filename,'r');
    while ~feof(fileID)
        [type, linedata] = readline(fileID);
        
        switch type
        case { 'blank', 'unknown'}
            continue;
        case 'layer'
            k = k + 1;
            laydata = linedata;
            if isfield(linedata, 'thickness')
                laydata.thickness = siParse(linedata.thickness);
            else
                j = j + 1;
                linedata.name = ['t' laydata.name];
                parameters(j) = parseParameter(linedata);
                laydata = rmfield(laydata, {'default', 'unit', 'plus', 'minus', 'plusminus'});
                laydata.thickness = parameters(j).value;
            end
            matid   = laydata.material;
            matdata = HBAR_materials(matid);
            config.layers(k)    = laydata;
            config.materials(k) = matdata;
        case 'contact'
            config.contact.resistance = siParse(linedata.resistance);
            config.contact.inductance = siParse(linedata.inductance);
        case 'radii'
            config.geometry.R1 = siParse(linedata.R1);
            config.geometry.R2 = siParse(linedata.R2);
            config.geometry.R3 = siParse(linedata.R3);
        case 'deltaR'
            config.geometry.deltaR = siParse(linedata.deltaR);
        case 'parameter'
            j = j + 1;
            parameters(j) = parseParameter(linedata);                    
        end
    end 
    fclose(fileID);
end

function [type, linedata] = readline(fileID)
    line = fgetl(fileID);
    %fprintf('%s\n', line);

    if isempty(line)
        type = 'blank';
        linedata = {};
        return;
    end
    
    % remove comments
    line = regexprep(line, '\s*#.*$', '');
    
    epmval = '(?<default>\d*\.?\d*)(?<unit>[pnumkMG]?\w*)\s*(?<plus>+\d*\.?\d*)?\s*(?<minus>-\d*\.?\d*)?\s*(?<plusminus>+-\d*\.?\d*)?';
    
    elayer = 'layer (?<type>[TPB][ELS]) (?<name>\w+) (?<material>\w+) (?<thickness>\d+[pnum]?m) "(?<description>[\w\s]+)";?';
    elayr1 = ['layer (?<type>[TPB][ELS]) (?<name>\w+) (?<material>\w+) ' epmval ' "(?<description>[\w\s]+)";?'];
    eradii = 'radii (?<R1>\d+um) (?<R2>\d+um) (?<R3>\d+um);?';
    edelta = 'deltaR (?<deltaR>d*\.?\d*um);?';    
    ecntct = 'contact (?<resistance>\d*\.?\d*[umkMG]?Ohm) (?<inductance>\d*\.?\d*[fpnu?]?H);?';
    eparam = ['parameter (?<name>\w+)\s+'  epmval];
    
    exprs = {elayr1 elayer eradii edelta ecntct eparam};
    exprs = strrep(exprs, ' ', '\s+');
    types = {'layer' 'layer' 'radii' 'deltaR' 'contact' 'parameter'};
    
    n = length(exprs);
    for k = 1:n
        linedata = regexp(line, exprs{k}, 'names');
        if ~isempty(linedata)
            type   = types{k};
            return;
        end
    end

    type = 'unknown';
    linedata = {};
end

function parameter = parseParameter(linedata)
    
    unit = linedata.unit;
    
    default = siParse([linedata.default unit]);
    plus = 0;
    minus = 0;
    
    if length(linedata.minus) > 1
        vstring = linedata.minus(2:end);
        minus = siParse([vstring unit]);
    end
    if length(linedata.plus) > 1
        vstring = linedata.plus(2:end);
        plus = siParse([vstring unit]);
    elseif length(linedata.plus) == 1
        plus = minus;
    end
    
    parameter.name = linedata.name;
    parameter.value = default;
    parameter.ubound = default + plus;
    parameter.lbound = default - minus;
end
