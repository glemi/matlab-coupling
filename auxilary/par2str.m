% converts a sequence of name-/value pairs into a string
% 
function str = par2str(args, pairsep, nvsep)
    opt pairsep char ';';
    opt nvsep char '=';
  
    n = length(args);
    names  = args(1:2:n-1);
    values = args(2:2:n);
    
    str = '';
    for k = 1:n/2
        name = names{k};
        value = strtrim(siPrefix(values{k}));
        attach = sprintf('%s%s%s%s', name, nvsep, value, pairsep);
        str = [str attach];
    end

    str = str(1:end-length(pairsep));
end

