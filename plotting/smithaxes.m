function smithaxes
    % Set XY data
    
    ax = gca;
    cla;
    
    gridlinepar = {'Color', [1 1 1]*.4, 'LineWidth', .5};
    
    [x, y] = mkgrid;
    angle = linspace(0, 2*pi, 129);
    hp = patch(sin(angle), cos(angle), 'w'); 
    hg = line(x, y, gridlinepar{:});
    hl = line([-1 1], [0 0], gridlinepar{:});
    
    skiplegend(hp, hg, hl);
    
    ax.XScale = 'lin';
    ax.YScale = 'lin';
    
    axis equal;
    axis off;
    
%     label(h); % rfchart.smith.label
end

function [xgrid, ygrid] = mkgrid()

    gridpoints =  [0.2000    0.5000    1.0000    2.0000    5.0000; ...
                   1.0000    2.0000    5.0000    5.0000   30.0000];

    % Initialize line data vectors
    N = 100; M = 128;
    nn = size(gridpoints,2);
    X = zeros(1,((2+4*nn)*(N+1))+M+1);
    Y = zeros(1,((2+4*nn)*(N+1))+M+1);

    % Full circle at R=50
    r = 50;
    x0 = r/(r+1);
    r0 = 1/(r+1);
    t = 0:M;
    s1 = 1;
    s2 = s1+M;
    X(1,s1:s2) = x0+r0*sin(t*2*pi/M);
    Y(1,s1:s2) = r0*cos(t*2*pi/M);

    % Full arcs at X = -50,+50
    s1 = s2+1;
    s2 = s1+2*N+1;
    z = z2g((linspace(0,50,N)).^2+1i*50*ones(1,N));
    X(1,s1:s2) = [NaN real(z) NaN real(conj(z))];
    Y(1,s1:s2) = [NaN imag(z) NaN imag(conj(z))];

    % R+jX grid
    for idx=1:nn
        s1 = s2+1;
        s2 = s1+4*N+3;
        r  = gridpoints(1,idx);
        x  = (linspace(0,sqrt(gridpoints(2,idx)),N)).^2;
        ZR = z2g(r*ones(1,N)+1i*x);
        ZX = z2g(x+1i*r*ones(1,N));
        X(1,s1:s2) = [NaN real(ZR) NaN real(conj(ZR)) NaN real(ZX)          ...
            NaN real(conj(ZX))];
        Y(1,s1:s2) = [NaN imag(ZR) NaN imag(conj(ZR)) NaN imag(ZX)          ...
            NaN imag(conj(ZX))];
    end

    % Record new data
    xgrid = X;
    ygrid = Y;
end


%%%%%%%%%%%%%%%%%%%
function g = z2g(z)
    % Z2G  Convert normalized impedance to reflection coefficient
    g = (z-1)./(z+1);
end
