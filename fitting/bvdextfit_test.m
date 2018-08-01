function fit = bvdextfit_test(f, Z, varargin)
    
    M = Z2M(Z);

    global debug_plot;
    if debug_plot
        progressPlot(f, M, true);
    end

    [~, fr] = findpeaks(phase(Z), f, 'NPeaks', 1, 'MinPeakDistance', diff(f([1 end]))/2);
    %fr = 2e9;
    Rm = 20;
    Qm = 4;
    C0 = abs(mean(1./(2*pi*f(end).*Z(end))));
    R0 = 1e-10; % not used
    Rp = 1e6; 
    Rc = real(Z(end));
    Lc = 0.1e-12;
    U  = 20e-12; % unknown parameter
    td = 0.04;

    initial  = [fr Rm Qm C0 R0 Rp Rc Lc U td];
    sinitial = v2c(fr, Rm, Qm, C0, R0, Rp, Rc, Lc, U, td);
    model    = mkmodel([fr Rm Qm C0 R0 Rp Rc Lc U td]);
    
    warning off MATLAB:rankDeficientMatrix;
    
    options  = statset('nlinfit');
    final    = nlinmatrixfit(f, M, model, sinitial, options, varargin{:});
    
    [fr, Rm, Qm, C0, R0, Rp, Rc, Lc, U, td] = c2v(initial, final);
    
    
    fit.Zfit  = M2Z(model(final, f));
    fit.final = c2v(initial, final);
    
    fit.ZfitModel = @(f)M2Z(model(final, f));
    
    final_clean = final;
    final_clean(5:8) = [0 1 0 1]*1e50;
    fit.final_clean = c2v(initial, final_clean);
    fit.Zfit_clean = M2Z(model(final_clean, f));
    
    fit.fr = fr;
    fit.Rm = Rm;
    fit.Qm = Qm;
    fit.C0 = C0;
    fit.R0 = R0;
    fit.Rp = Rp;
    fit.Rc = Rc;
    fit.Lc = -Lc;    
    fit.U  = U;
    fit.td = td;
    
    fit.Cm = 1/(2*pi*fr*Qm*Rm);
    fit.Lm = 1/(2*pi*fr)*Qm*Rm;
    fit.fa = sqrt(1/(fit.Lm*fit.Cm)*(1 + fit.Cm/C0))/(2*pi);
end

function M = Z2M(Z)
    Z = Z(:);

%     M(:,1) = real(Z);
%     M(:,2) = log(-imag(Z));
    M(:,1) = abs(Z);
    M(:,2) = phase(Z);
end

function Z = M2Z(M)    
    [re, im] = pol2cart(M(:,2), M(:,1));
    Z = complex(re, im);
    %Z = complex(M(:,1), -exp(M(:,2)));
end


function Z = bvdmodel(f, fr, Rm, Qm, C0, R0, Rp, Rc, Lc, U, td)
    Cm = 1/(2*pi*fr*Qm*Rm);
    Lm = (Qm*Rm)/(2*pi*fr);
    
    ZC0 = 1./(2i*pi*f*C0*exp(-1i*td)); 
    ZCm = 1./(2i*pi*f*Cm);
    ZLm =     2i*pi*f*Lm;
    ZLc =     2i*pi*f*Lc;
    
    par = @(x,y)x.*y./(x + y);
    
    Zm = ZLm + ZCm + Rm;
    Z0 = ZC0; %+ R0;
    Z = par(Zm, Z0);
    %Z = Z + Rc - ZLc + 1./(U*2*pi*f);
    Z = par(Rp, Z) + Rc - ZLc; %+ 1./(U*2*pi*f);
    %Z = Z + Rc - ZLc;  %+ 1./(U*2*pi*f).^2;
end

function model = mkmodel(initial)
    model = @(c,f)wrapper(c,f);
    
    function M = wrapper(coeffs, f)
        [fr, Rm, Qm, C0, R0, Rp, Rc, Lc, U, td] = c2v(initial, coeffs);    
        Z = bvdmodel(f, fr, Rm, Qm, C0, R0, Rp, Rc, Lc, U, td);
        M = Z2M(Z);
        global debug_plot;
        if debug_plot
            progressPlot(f, M, false);
        end    
    end
end

function progressPlot(f, M, keep)
    persistent h1 h2; try delete([h1 h2]); end; %#ok
    fig bvdextfit;
    if keep
        clf;
        subplot(2,1,1);      plot(f, M(:,1),  '-', 'LineWidth', 0.5);
        subplot(2,1,2);      plot(f, M(:,2),  '-', 'LineWidth', 0.5);
    else
        subplot(2,1,1); h1 = plot(f, M(:,1), 'k-', 'LineWidth', 0.5);
        subplot(2,1,2); h2 = plot(f, M(:,2), 'k-', 'LineWidth', 0.5);
    end
    drawnow;
end


function coeffs = v2c(fr, Rm, Qm, C0, R0, Rp, Rc, Lc, U, td)
    coeffs = [fr Rm Qm C0 R0 Rp Rc Lc U td];
    
    scalefact = 10.^(-floor(log10(coeffs)));
    coeffs = abs(scalefact.*coeffs);
end

function [fr, Rm, Qm, C0, R0, Rp, Rc, Lc, U, td] = c2v(initial, coeffs)
    scalefact = 10.^(floor(log10(initial)));
    coeffs = abs(scalefact.*coeffs);
    
    fr = coeffs(1);
    Rm = coeffs(2); 
    Qm = coeffs(3);
    C0 = coeffs(4);
    R0 = coeffs(5);
    Rp = coeffs(6);
    Rc = coeffs(7);
    Lc = coeffs(8);
    U  = coeffs(9);
    td = coeffs(10);
end

