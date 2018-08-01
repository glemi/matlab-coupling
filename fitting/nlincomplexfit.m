function [final,R,J,COVB,MSE] = nlincomplexfit(x, y, model, initial, varargin)
    xx = [x(:); x(:)];
    yy = [real(y(:)); imag(y(:))];
        
    function yy = wrapper(coeffs, ~)
        Y = model(coeffs, x);
        yy = [real(Y(:)); imag(Y(:))];
    end

    i = find(strcmp('Weights', varargin), 1);
    if ~isempty(i)
        w = varargin{i+1};
        varargin{i+1} = [w(:); w(:)];
    end
    
    [final,R,J,COVB,MSE] = nlinfit(xx, yy, @wrapper, initial, varargin{:});
end

