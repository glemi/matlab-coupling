function [Zfit, data] = hbarmultifit_v3(f, z, config, cnames, curves, params, varargin)
        
    if ischar(config) && exist(config, 'file')
        config = HBAR_loadconfig(config);
    elseif ~isstruct(config)
        error 'Invalid HBAR config parameter supplied';
    end
    
    if nargin < 4
        cnames =  {'cPiezo' 'ePiezo' 'epsPiezo' 'tPiezo' 'tBotEl' 'tTopEl' };
    end
    if nargin < 5
        params = {};
    end

    options = {'MaxIter', 20, 'FunValCheck', 'off', 'DerivStep', 0.1, ...
        'RobustWgtFun', 'bisquare'};
    options = statset(options{:});
    
    f = [f(:)];
    N = 400;
    [M, F] = HBAR_postprocess(f, z, curves, N);  
    fig mplot; clf; HBAR_Mplot(F, M, curves);
    varargin = adjustweights(varargin,N);
    
    initial = HBAR_getconfigvalues(config, cnames);
    sinitial = scale(initial, initial);
    model = mkmodel(f, cnames, config, curves, params, initial, N);
    
    [sfinal, R,J,C,E] = nlinmatrixfit(model, sinitial, F, M, options);
    [Mfit, Zfit, hbardata] = model(sfinal, f);
    
    final  = unscale(initial, abs(sfinal));
    
    sj = sum(J,1)/length(J);
    
    cpairs = nvpairs(cnames, final);
    config = HBAR_parameters(config, cpairs);
    
    data.f       = f;
    data.F       = F;
    data.Zmeas   = z;
    data.Zfit    = Zfit;
    data.Mfit    = Mfit;
    data.config  = config;
    data.curves  = curves;
    data.names   = cnames;
    data.initial = initial;
    data.final   = final;
    data.hbar    = hbardata;
    data.model   = @(f)model(final, f);
    
    data.Residuals = R;
    data.Jacobian = J;
    data.Covariance = C;
    data.ErrorVariance = E;
end

function args = adjustweights(args, N)
    i = find(strcmp('Weights', args), 1);
    if ~isempty(i)
        w = args{i+1};
        n = length(w);
        s = round(n/N);
        args{i+1} = w(1:s:n);
    end
end

function model = mkmodel(f, cnames, config, spectra, params, initial, N)
    opt N double length(f);

    function [M, Z, data] = hbarmodel(coeffs, ~)
        persistent lastcoeffs lastZ lastM lastData h;
        if all(size(lastcoeffs) == size(coeffs)) && all(coeffs == lastcoeffs)
            fprintf('same coefficients!\n');
            M = lastM;
            Z = lastZ;
            data = lastData;
        else
            %fprintf('%.4f ', coeffs); fprintf('\n');
            coeffs = unscale(initial, coeffs);
            cpairs = nvpairs(cnames, abs(coeffs));

            % in this order, coeffs override fixed parameters:
            [Z, data] = HBAR_v3(f, config, params{:}, cpairs{:});
            [M, F]    = HBAR_postprocess(f, Z, spectra, N);

            try delete(h); end %#ok
            fig mplot; 
            h = HBAR_Mplot(F, M, spectra, 'k-');
            drawnow;
            
            lastcoeffs = coeffs; lastZ = Z; lastM = M; lastData = data; 
        end     
    end
    
    model = @hbarmodel;
end

function scaled = scale(initial, coeffs)
    scalefact = 10.^(-floor(log10(initial)));
    scaled = scalefact.*coeffs;
end

function unscaled = unscale(initial, scaled)
    scalefact = 10.^(floor(log10(initial)));
    unscaled = scalefact.*scaled;
end

