function assignargs(varargin)
    
    n = length(varargin);
    if mod(n,2) > 0
        error 'need even number of arguments as name/value pairs';
    end
    
    for k = 1:2:n
        name = varargin{k};
        value = varargin{k+1};
        
        if isvarname(name)
            assignin('caller', name, value);
        end
    end
end