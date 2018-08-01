function [final,R,J,COVB,MSE] = nlinmultifit(X, Y, models, initial, varargin)

	[m, n] = size(Y);
	[Y, s] = scale(Y);
	
    xx = repmat(X, n, 1);
    yy = reshape(Y, n, 1);
	
	model = models(1);
	
    function yy = wrapper(coeffs, xx)
		yy = zeros(size(xx));
		for k = 1:n
			range = [1:m] + (k-1)*m;
			model = models(k);
			x = xx(range);
			y = model(coeffs, x);
			yy(range) = scale(y, s(:,k));
		end
    end
    
	[final,R,J,COVB,MSE] = nlinfit(xx, yy, @wrapper, initial, varargin{:});
end

function [y, s] = scale(y, s)
	m = size(y, 1);
	
	if nargin >= 2
		a = s(1,:);
		b = s(2,:);
	else
		miny = min(y, 1);
		maxy = max(y, 1);
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