% function pairs = nvpairs(names, values)
% 
% Creates a cell array of name/value pairs that can be used as input for
% various functions, such as plot, figure, etc. 
% 
% usage:
%           names = {'LineWidth' 'MarkerSize'};
%           values = [2, 10];
%           pairs = nvpairs(names, values);
%           plot(x, y, pairs{:});
%           
%           names = {'LineStyle' 'LineWidth' 'MarkerSize'};
%           values = {'.-' 2 10};
%           pairs = nvpairs(names, values);
%           plot(x, y, pairs{:});
%           
function pairs = nvpairs(names, values)
    
    if isnumeric(values)
        values = num2cell(values);
    end

    pairs = [names(:) values(:)]';
    pairs = pairs(:)';    
end

