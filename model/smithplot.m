function hline = smithplot(varargin)
    hline = plot(varargin{:});
    
    x = hline.XData;
    y = hline.YData;
    z = complex(x, y);
    
    g = z2gamma(z, 50);
    x = real(g);
    y = imag(g);
    
    hline.XData = x;
    hline.YData = y;
    
    drawnow;
end

