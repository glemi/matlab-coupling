function [zfit, data] = hbarmultifit(f, z, config, cnames, params, varargin)
        
    if ischar(config) && exist(config, 'file')
        config = HBAR_loadconfig(config);
    elseif ~isstruct(config)
        error 'Invalid HBAR config parameter supplied';
    end
    
    if nargin < 4
        cnames =  {'cPiezo' 'ePiezo' 'epsPiezo'   'tPiezo' 'tBotEl' 'tTopEl' };
        initial = [380e9     1.5       12          880e-9   250e-9   150e-9   ];
    else
        initial = HBAR_getconfigvalues(config, cnames);
    end
    
    if nargin < 5
        params = {};
    end
    
    options = {'MaxIter', 15, 'FunValCheck', 'off', 'DerivStep', 0.005, ...
        'RobustWgtFun', 'bisquare'};
    
    f = [f(:)];
    N = 200;
    [M, F] = z2M(f, z, N);  
    fig mplot; clf; Mplot(F, M, f, z, 'first');
    
    sinitial = scale(initial, initial);
    fine   = mkmodel(f, cnames, config, params, initial);
    coarse = mkmodel(f, cnames, config, params, initial, N);
    
    [sfinal, res, jac] = nlinmatrixfit(F, M, coarse, sinitial, options{:}, varargin{:});
    sfinal = abs(sfinal);
    sj = sum(jac,1)/length(jac);
    
    [M, Z, hbardata] = fine(sfinal, f);
    
    zfit = Z;
    
    data.names   = cnames;
    %data.units   = cunits;
    data.initial = initial;
    data.final   = unscale(initial, sfinal);
    data.hbar    = hbardata;
    data.model   = @(f)coarse(final, f);
end

function model = mkmodel(f, cnames, config, params, initial, N)
    opt N double length(f);

    function [M, Z, data] = hbarmodel(coeffs, ~)
        persistent lastcoeffs;
        if all(size(lastcoeffs) == size(coeffs)) && all(coeffs == lastcoeffs)
            fprintf('same coefficients!\n');
        else
            lastcoeffs = coeffs;
        end
        
        fprintf('%.4f ', coeffs); fprintf('\n');
        coeffs = unscale(initial, coeffs);
        values = num2cell(abs(coeffs));
        cpairs = [cnames; values];
        
        [Z, data] = HBAR_v3(f, config, cpairs{:}, params{:});
        [M, F] = z2M(f, Z, N);
        [~, ru, rl] = deripple(f, real(Z), N);
        [~, iu, il] = deripple(f, imag(Z), N);
        Zu = complex(ru, iu);
        Zl = complex(rl, il);
        Mplot(F, M, F, [Zu Zl]);
    end
    
    model = @hbarmodel;
end

function [M, F] = z2M(f, Z, N)
    if nargin < 3
        N = [];
    end
    
    rippledata = multiRipplefit(f, abs(Z), 3);
    fr   = [rippledata.fr];
    dfr  = diff([rippledata.fr]);
    keff = [rippledata.kt2];
    C0   = [rippledata.C0];
    Cm   = [rippledata.Cm];
    %Rm   = [rippledata.R];
    
    F = linspace(f(1), f(end), N)';
    dfr  = interp1(fr(2:end), dfr, F);
    keff = interp1(fr, keff, F);
    C0   = interp1(fr, C0, F);
    Cm   = interp1(fr, Cm, F);
    %Rm   = interp1(fr, Rm, F);
    
    M = [dfr keff Cm];
    if any(any(~isfinite(M)))
        warning('model is generating non-finite numbers');
    end
end

function scaled = scale(initial, coeffs)
    scalefact = 10.^(-floor(log10(initial)));
    scaled = scalefact.*coeffs;
end

function unscaled = unscale(initial, scaled)
    scalefact = 10.^(floor(log10(initial)));
    unscaled = scalefact.*scaled;
end

function Mplot(F, M, f, Z, options)
    persistent p;
    style = 'k-';
    
    fig mplot;
    if nargin > 4 && strcmp(options, 'first');
        subplot(3,2,1); title dfr;   plot(F, M(:,1));
        subplot(3,2,2); title keff;  plot(F, M(:,2));
        subplot(3,2,3); title Cm;    plot(F, M(:,3));
        %subplot(3,2,4); title Rm;    plot(F, M(:,4));
        subplot(3,2,5); title |Z|;   plot(f, abs(Z)); 
        subplot(3,2,6); title \phi;  plot(f, phase(Z)); 
        yscale log; xscale log;
    else 
        try delete(p); end;
        subplot(3,2,1); p(1,:) = plot(F, M(:,1), style);
        subplot(3,2,2); p(2,:) = plot(F, M(:,2), style);
        subplot(3,2,3); p(3,:) = plot(F, M(:,3), style);
        %subplot(3,2,4); p(4,:) = plot(F, M(:,4), style);
        subplot(3,2,5); p(5,:) = plot(f, abs(Z(:,1)), style);
        subplot(3,2,5); p(6,:) = plot(f, abs(Z(:,2)), style);
        subplot(3,2,6); p(7,:) = plot(f, phase(Z(:,1)), style);
        subplot(3,2,6); p(8,:) = plot(f, phase(Z(:,2)), style);
    end
    drawnow;
end
