% function struct = var2struct( varargin )
%
% Makes all variables passed as arguments a field in a structure, where the
% fieldnames correspond to the variable names. Uses inputname to determine
% the variable name of the arguments. 
%
% usage: 
%           
%
%
% see also: inputname dissolve
function struct = var2struct( varargin )
    n = nargin;
    
    for k = 1:n
        name = inputname(k);
        if isvarname(name)
            struct.(name) = varargin{k};
        end
    end
end

