classdef HBAR_Fit < handle
    % class HBAR_Fit 
    %
    % usage (example):
    %           f = data.f;
    %           Zmeas = squeeze(data.z(1,1,:));
    %           Zclean = HBAR_removeparasitics_v2(f, Zmeas);
    % 
    %           Hfit = HBAR_Fit;
    %           Hfit.FitOptions.MaxIterations = 10;
    %           Hfit.FitOptions.FiniteDifferenceStepSize = 0.01;
    %           Hfit.FitParams = {'cPiezo' 'ePiezo' 'epsPiezo' 'tPiezo' 'tBotEl' 'tTopEl' 'QSubst'};
    %           Hfit.ParamErrs = [20e9       0.5       5         40e-9   20e-9    60e-9    500];
    %           Hfit.Curves = {'ripples:keff' 'ripples:Cm' 'ripples:diff(fr)'};;
    %           Hfit.execute(f, Zclean);
    
    properties
        Config = HBAR_loadconfig('HBAR_config.txt');
        Curves = {'ripples:keff' 'ripples:Cm' 'ripples:diff(fr)'};
        FitParams = {'cPiezo' 'ePiezo' 'epsPiezo' 'tPiezo' 'tBotEl' 'tTopEl' };
        ParamErrs = [20e9 0.5 1 20e-9 20e-9 60e-9];
        FixedParams = {};
        Npoints = 300;
        
        FitOptions = optimoptions('lsqcurvefit');
        SmoothDfr = true;
        PlotEnable = true;
    end
    
    properties (SetAccess = private, GetAccess = public)
        f
        F
        Zmeas
        Mmeas
        Zfit 
        Mfit 
        HBARfit
        InitialCoeffs 
        InitialConfig     % includes fixed parameters (c33E, epsf, etc)
        InitialProcConfig % preprocessed initl config (c33D, epsS, kt2 etc)
        FinalCoeffs
        FinalConfig       % final config (c33D, epsS, kt2 etc)
        FinalProcConfig   % preprocessed final config (c33D, epsS, kt2 etc)
        FitOutput
        
        fitdata % compatibility with hbarmultifit_v2;
    end
    
    methods
        function this = HBAR_Fit(config)
            if nargin > 0 
                if isstruct(config)
                    this.Config = config;
                elseif ischar(config) && logical(exist(config, 'file'))
                    this.Config = HBAR_loadconfig(config);
                else
                    error 'Invalid config parameter!';
                end
            end
        end
        
        function execute(this, f, Z)
            f = f(:);
            
            N      = this.Npoints;
            config = this.Config;
            cnames = this.FitParams;
            params = this.FixedParams;
            curves = this.Curves;
            errors = this.ParamErrs;
            
            if this.SmoothDfr
                curves1 = strrep(curves, 'ripples:dfr', 'ripples:smooth(dfr)');
                curves1 = strrep(curves1, 'ripples:diff(fr)', 'ripples:smooth(dfr)');
            else
                curves1 = curves;
            end
            
            iconfig = HBAR_parameters(config, params);
            ipconfig = proc_config(iconfig);
            
            [M, F]  = HBAR_postprocess(f, Z, curves1, N, true);  
            initial = HBAR_parameters(iconfig, cnames);
            model   = mkmodel(f, cnames, iconfig, curves, params, initial, N, this.PlotEnable);
            
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
            fconfig = HBAR_parameters(iconfig, cpairs);
            % replace c33E -> c33D, epsf -> epsS, (kt,c,eps) ~~> e33 
            fpconfig = proc_config(fconfig);
            ipconfig = proc_config(iconfig);
            
            this.f     = f;
            this.F     = F;
            this.Zmeas = Z;
            this.Mmeas = M;
            this.Zfit  = Zfit; %#ok
            this.Mfit  = Mfit; %#ok
            this.HBARfit = hbardata;
            this.InitialCoeffs = initial;
            this.InitialConfig = iconfig;
            this.InitialProcConfig = ipconfig;
            this.FinalCoeffs = final;
            this.FinalConfig = fconfig;
            this.FinalProcConfig = fpconfig; % preprocessed config
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
            
            data.config_initial = this.Config;
            data.config_final   = this.FinalConfig;

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
            [Z, data] = HBAR_v4(f, config, params{:}, cpairs{:});
            [M, F]    = HBAR_postprocess(f, Z, spectra, N);

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

function pconfig = proc_config(config)
    pconfig = config;
    hbar = HBAR_preprocess_v2(config);
    layers = [hbar.layers];
    pconfig.materials = [layers.material];
end

function scaled = scale(initial, coeffs)
    scalefact = 10.^(-floor(log10(initial)));
    scalefact(coeffs == 0) = 1;
    scaled = scalefact.*coeffs;
end
function unscaled = unscale(initial, scaled)
    scalefact = 10.^(floor(log10(initial)));
    scalefact(initial == 0) = 1;
    unscaled = scalefact.*scaled;
end