classdef HBAR_HFit < handle
    % Class HBAR_HFit
    % 
    % usage: 
    % 
    % 
    % 
    
    properties
        Config = HBAR_loadconfig('HBAR_config.txt');
        Curves = {'ripples:keff' 'ripples:Cm' 'ripples:diff(fr)'};
        FitParams = {'cPiezo' 'ePiezo' 'epsPiezo' 'tPiezo' 'tBotEl' 'tTopEl' };
        ParamErrs = [20e9 0.5 1 20e-9 20e-9 60e-9];
        FixedParams = {};
        Npoints = 300;
        
        FitOptions = optimoptions('lsqcurvefit');
        
        PlotEnable = true;
    end
    
    properties (SetAccess = private, GetAccess = public)
        f
        F
        Zmeas
        Mmeas
        Hmeas
        Zfit 
        Mfit 
        Hfit
        HBARfit
        InitialCoeffs
        FinalCoeffs
        FinalConfig
        FitOutput
        
        fitdata % compatibility with hbarmultifit_v2;
    end
    
    methods
        function execute(this, f, Z)
            f = f(:);
            
            N      = this.Npoints;
            config = this.Config;
            cnames = this.FitParams;
            params = this.FixedParams;
            curves = this.Curves;
            errors = this.ParamErrs;
            
            bvd = bvdfit(f,Z);
            H = (2i*pi*f*bvd.C0).*Z;
            
            [M, F]  = HBAR_postprocess(f, H, curves, N, true);  
            initial = HBAR_getconfigvalues(config, cnames);            
            model   = mkmodel(f, cnames, config, curves, params, initial, N, this.PlotEnable);
            
            if this.PlotEnable
                fig mplot; clf; 
                HBAR_Mplot(F, M, curves);
            end
            
            sinitial = scale(initial, initial);
            options  = this.FitOptions;
            
            if length(errors) == length(cnames)
                slb = scale(initial, initial - errors);
                sub = scale(initial, initial + errors);
                [sfinal, output] = lsqmatrixfit(model, sinitial, F, M, slb, sub, options);
            else
                [sfinal, output] = lsqmatrixfit(model, sinitial, F, M, options);
            end
            
            
            [Mfit, Zfit, hbardata] = model(sfinal, f);  %#ok
            final  = unscale(initial, abs(sfinal));

            cpairs = nvpairs(cnames, final);
            config = HBAR_parameters(config, cpairs);
            
            this.f     = f;
            this.F     = F;
            this.Zmeas = Z;
            this.Mmeas = M;
            this.Hmeas = H;
            this.Zfit  = Zfit; %#ok
            this.Hfit  = hbardata.H;
            this.Mfit  = Mfit; %#ok
            this.HBARfit = hbardata;
            this.InitialCoeffs = initial;
            this.FinalCoeffs = final;
            this.FinalConfig = config;
            this.FitOutput = output;
            this.update_fitdata;
        end
        
    end
    
    methods (Access = private)
        function update_fitdata(this)
            % for compatibility with hbarmultifit_v2;
            data.f       = this.f;
            data.F       = this.F;
            data.Zmeas   = this.Zmeas;
            data.Zfit    = this.Zfit;
            data.Mfit    = this.Mfit;
            data.config  = this.FinalConfig;
            data.curves  = this.Curves;
            data.names   = this.FitParams;
            data.initial = this.InitialCoeffs;
            data.final   = this.FinalCoeffs;
            data.hbar    = this.HBARfit;
            data.model   = @(f)model(final, f);

            data.Residuals   = this.FitOutput.residual;
            data.Jacobian    = this.FitOutput.jacobian;
            data.Covariance  = [];
            data.ErrorVariance = [];
            this.fitdata = data;
        end

    end
    
end

function model = mkmodel(f, cnames, config, spectra, params, initial, N, plotEnable)
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
            [M, F]    = HBAR_postprocess(f, data.H, spectra, N);

            if plotEnable
                try delete(h); end %#ok
                fig mplot; 
                h = HBAR_Mplot(F, M, spectra, 'k-');
                drawnow;
            end

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