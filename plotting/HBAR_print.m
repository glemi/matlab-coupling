% outtext = HBAR_print(mode, varargin)
% 
% Prints selected parameters of a HBAR configuration in human-readable
% format, either as plain text or latex encoded. 
%
% usage: 
%           % output into command window
%           % name-value pairs as vararing
%           HBAR_print('plain', 'name1', value1, 'name2', 'value2'...);
%
%           % names as cellstring and values as array
%           HBAR_print('plain', {'name1' 'name2' ...}, [value1 value2 ...]);
%
%           % names as cellstring and config structure
%           HBAR_print('plain', {'name1' 'name2' ...}, config);
%           
%           % output to variable as plain text
%           str = HBAR_print('plain', ...);
%           
%           % output to variable as latex encoded text
%           str = HBAR_print('latex', ...);
%
%           % output to variable as latex encoded 
%           str = HBAR_print('blocklatex', ...);
%
% see also:
%           HBAR_parameters, HBAR_loadconfig        
%   
function outtext = HBAR_print(mode, varargin)
    [codes, values] = parseargs(varargin{:});
    
    n = length(codes);
    outtext = {};
    
    switch mode
        case 'plain'
            for k = 1:n
                partext = printpar(codes{k}, values(k), mode);
                outtext{k} = partext;
            end
        case 'latex'
            for k = 1:n
                partext = printpar(codes{k}, values(k), mode);
                outtext{k} = ['$' partext '$'];
            end
        case 'blocklatex'
            outtext = '$\begin{array}{rl}';
            for k = 1:n
                partext = printpar(codes{k}, values(k), mode);
                outtext = [outtext partext];
            end
            outtext = [outtext ' \end{array}$'];
    end
    
    if nargout == 0
        outtext = cellstr(outtext);
        outtext = strjoin(outtext, '\n');
        disp(outtext);
    end
    
end

function partext = printpar(code, value, mode)
    
    config = getprintconfig(code);
    
    switch mode
        case 'plain'
            if isempty(config.unit)
                valuetext = sprintf(config.format, value);
            elseif strcmp(config.unit, '%')
                valuetext = sprintf([config.format '%%'], value*100);
            else
                valuetext = siPrefix(value, config.unit, config.format);
            end
            partext = sprintf('%8s: %4s', config.symbol, valuetext);
            
        case 'latex'
            if isempty(config.unit)
                valuetext = sprintf(config.format, value);
            elseif strcmp(config.unit, '%')
                valuetext = sprintf([config.format '\\%%'], value*100);
            else
                valuetext = siPrefix(value, config.lunit, config.format);
                valuetext = regexprep(valuetext, '(\d*\.?\d*)\s?(\D+)', '$1\\,$2');
            end
            partext = sprintf('%8s: \\rm %4s', config.lsymbol, valuetext);
            partext = strrep(partext, char(181), '\mu{}');
            partext = strrep(partext, char(65533), '\mu{}');
            
        case 'blocklatex'
            if isempty(config.unit)
                valuetext = sprintf(config.format, value);
            elseif strcmp(config.unit, '%')
                valuetext = sprintf([config.format '\\%%'], value*100);
            else
                valuetext = siPrefix(value, config.lunit, config.format);
                valuetext = regexprep(valuetext, '(\d*\.?\d*)\s?(\D+)', '$1\\,$2');
            end
            partext = sprintf('%8s:& \\rm %4s \\\\', config.lsymbol, valuetext);
            partext = strrep(partext, char(181), '\mu{}');
            partext = strrep(partext, char(65533), '\mu{}');
    end
end

function pconfig = getprintconfig(code)
    parcodes = {'eps' 'td' 'nu' 'd' 'Q' 'c' 'e' 'r' 't'};
    laycode = '';
    parcode = code;
    
    n = length(parcodes);
    for k = 1:n
        pc = parcodes{k};
        m = length(pc);
        if strncmp(code, pc, m)
            laycode = code(m+1:end);
            parcode = pc;
            break;
        end
    end
    
    s_eps = '\varepsilon_{33}'; 
    u_r   = '\Omega\cdot{m}';

    switch parcode
		case 'R1', lsymbol='R_1';     unit='m';     lunit='m';      format='%.0f';
		case 'R2', lsymbol='R_2';     unit='m';     lunit='m';      format='%.0f';
		case 'R3', lsymbol='R_3';     unit='m';     lunit='m';      format='%.0f';
		case 'Rc', lsymbol='R_c';     unit='Ohm';   lunit='\Omega'; format='%3.*f';
		case 'Lc', lsymbol='L_c';     unit='H';     lunit='H';      format='%3.*f';
        case 'xC0',lsymbol='x_{C_0}'; unit='';      lunit='';       format='%.2f'; 
        case 'td', lsymbol='\tan(\delta)'; unit='%'; lunit='%';     format='%.2f';
		case 't',  lsymbol='t';       unit='m';     lunit='m';      format='%.0f';
		case 'd',  lsymbol='\rho';    unit='kg/m3'; lunit='kg/m^3'; format='%3.*f';
		case 'Q',  lsymbol='Q';       unit='';      lunit='';       format='%.0f';
		case 'nu', lsymbol='\nu';     unit='';      lunit='';       format='%.2f';
		case 'c',  lsymbol='c_{33}';  unit='Pa';    lunit='Pa';     format='%.0f';
		case 'e',  lsymbol='e_{33}';  unit='';      lunit='';       format='%.2f';
		case 'eps',lsymbol=s_eps;     unit='';      lunit='';       format='%.1f';
		case 'r',  lsymbol='\rho';    unit='Ohm?m'; lunit=u_r;      format='%3.*f';
        case 'C0', lsymbol='C_0';     unit='F';     lunit='F';      format='%3.*f';
        case 'kt2',lsymbol='k_t^2';   unit='%';     lunit='%';      format='%.2f';
        otherwise, error 'unknown parameter code';
    end    
    
    symbol = code;
    lsymbol = [lsymbol sprintf('\\,_{\\rm{}%s}', laycode)];
    lunit = [lunit];
    pconfig = var2struct(symbol, lsymbol, unit, lunit, format);
end

function [codes, values] = parseargs(varargin)
    args = varargin;
    n = nargin;
    
    if n == 2 && ischar(args{1}) && isnumeric(args{2})
        values = args{2};
        codes = repmat(args(1), size(values));
    elseif iscellstr(args(1:2:n)) && isnumeric([args{2:2:n}])
        % name-values pairs
        codes = args(1:2:n);
        values = [args{2:2:n}];
    elseif iscellstr(args{1}) && isnumeric(args{2}) && (length(args{1}) == length(args{2}))
        % cell array of codes in arg{1} and values in arg{2}
        codes = args{1};
        values = args{2};
    elseif iscellstr(args{1}) && isstruct(args{2})
        codes = args{1};
        config = args{2};
        values = HBAR_getconfigvalues(config, codes);
    elseif ischar(args{1}) && isstruct(args{2})
        codes = args{1};
        config = args{2};
        values = HBAR_getconfigvalues(config, codes);
    else
        error 'invalid input';
    end
end