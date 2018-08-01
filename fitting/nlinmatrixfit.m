%function [final,R,J,COVB,MSE] = nlinmatrixfit(X, Y, model, initial, opts, varargin)
function [final,R,J,COVB,MSE] = nlinmatrixfit(X, Y, model, initial, opts, varargin)
    opt opts struct statset('nlinfit');

	[m, n] = size(Y);
	[Y, s] = scale(Y);
	
    xx = repmat(X, n, 1);
    yy = reshape(Y, n*m, 1);
	
    i = find(strcmp('Weights', varargin), 1);
    if ~isempty(i)
        w = varargin{i+1};
        varargin{i+1} = repmat(w,n,1);
    end
    
	function yy = wrapper(coeffs, xx)
		x  = xx(1:m);
		Y  = model(coeffs, x);
		Y  = scale(Y, s);
		yy = reshape(Y, n*m, 1);
    end
    
	[final,R,J,COVB,MSE] = nlinfit(xx, yy, @wrapper, initial, opts, varargin{:});
    %[final,R,~,~,~,~,J] = lsqcurvefit(@wrapper, initial, xx, yy)
end

function [y, s] = scale(y, s)
	m = size(y, 1);
	
	if nargin >= 2
		a = s(1,:);
		b = s(2,:);
	else
		miny = min(y, [], 1);
		maxy = max(y, [], 1);
		a = maxy-miny;
		b = miny;
	end
		
	A = repmat(a, m, 1);
	B = repmat(b, m, 1);
		
	y = (y-B)./A;
	s = [a; b];
end

function y = unscale(y, s)
	a = s(1,:);
	b = s(2,:);
	
	A = repmat(a, m, 1);
	B = repmat(b, m, 1);
	
	y = A.*y + B;
end