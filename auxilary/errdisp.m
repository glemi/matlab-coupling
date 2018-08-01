function compact = errdisp(err, truncate)
    if nargin < 2
        caller = dbstack;
        truncate = caller(2).name;
    end

    s = getReport(err, 'extended');
    
    if ~isempty(err.identifier)
        compact = sprintf('\n%s\n(%s)\n', err.message, err.identifier);
    else
        compact = sprintf('\n%s\n\n', err.message);
    end
    
    lines = strsplit(s, '\n');
    n = length(lines);
    for k = 1:n
        line = lines{k};
        
        if ~any(strfind(line, 'Error in')) && ~any(strfind(line, 'Error using'))
            continue;
        end
        
        line = strrep(line, 'Error in ', '');
        line = strrep(line, 'Error using ', ''); 
        line = regexprep(line, '([^<])/', '$1.');
        line = strrep(line, '(varargin)', '');
        line = strrep(line, '(varargin{:})', '');
    	compact = [compact ' ' line 10]; %#ok;
        
        if any(strfind(line, truncate))
            break;
        end
    end
    
    if ~isempty(err.cause)
        if nargin == 2
            compact = [compact 10 'caused by:' errdisp(err.cause{1}, truncate)];
        else
            compact = [compact 'caused by' 10 errdisp(err.cause{1})];
        end
    end
    
    if nargout == 0
        fprintf(2, '%s\n', compact);
    end
end

