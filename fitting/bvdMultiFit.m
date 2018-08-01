function [data, zfit] = bvdMultiFit(f, z, n, varargin)
    hf = gcf;
    global debug_plot;
    
    i = strcmpi('recycle', varargin);
    if any(i)
        varargin = varargin(~i);
        initial = recycle([], f, z);
    else
        initial = coeff_estimate(f, z, n); 
    end
    opts    = statset(varargin{:});
    logf    = log(f);
    final   = nan(size(initial));   
 
    if ~any(isnan(initial))
        if debug_plot
            fig bvdfit:fitdebug; clf;
            plot(logf, log(abs(z)), 'LineWidth', 5);
            ylim(ylim.*[0.75 1.25]);
            axis manual;

            fig bvdfit:coeffs; clf;
            prepareCoeffPlot;

            final = nlinfit(logf, log(abs(z)), @(c,f)wrapper(n,c,f), initial, opts);
            logz = wrapper(n, final, logf);

            fig bvdfit:fitdebug;
            plot(logf, abs(logz), '--', 'LineWidth', 2);        
        else
            final = nlinfit(logf, log(abs(z)), @(c,f)wrapper(n,c,f), initial, opts);
            logz = wrapper(n, final, logf);
        end
        
        zfit = exp(logz);
        data = fill_data_struct(initial, final, n, f, zfit, z);
        %if data.valid
            recycle(final);
        %end
    else
        zfit = nan(size(f));
        data = fill_data_struct(initial, final, n, f, zfit, z);
        %if data.valid
            recycle(final);
        %end
    end
    
    
    figure(hf);
end

function ok = checkfr(initial, final)
    dfri = initial(2);
    dfrf = final(2);
    
    ok = abs(dfrf/dfri - 1) < 2e-3;
end

function outcoeffs = recycle(incoeffs, f, z)
    persistent coeffs;
    if nargin == 1
        coeffs = incoeffs;
    end
    if nargout == 1 && ~isempty(coeffs)
        [fr, dfr, Qm, Rm, C0, Rs] = c2v(coeffs);
        
        fcenter = mean(f([1 end]));
        njumps = round((fcenter - fr)/dfr);
        frnew = fr + njumps*dfr;
        outcoeffs = v2c(frnew, dfr, Qm, Rm, C0, Rs);
    end
end

function coeffs = v2c(fr, dfr, Qm, Rm, C0, Rs)
    fr  =  fr*1e-9;
    dfr = dfr*1e-6;
    Qm  =  Qm*1e-3;
    Rm  =  Rm*1e-1;
    C0  =  C0*1e12;
    coeffs = abs([fr dfr Qm Rm C0 Rs]);
end

function [fr, dfr, Qm, Rm, C0, Rs] = c2v(coeffs)
    n = (length(coeffs) - 4);

    fr  =     coeffs(1)*1e9;
    dfr =     coeffs(2)*1e6;
    Qm  =     coeffs(3)*1e3;
    Rm  =     coeffs(4)*1e1; 
    
    C0  =     coeffs(5)*1e-12;
    Rs  = abs(coeffs(6));
end

function logZ = wrapper(n, coeffs, logf)
    global debug_plot;
    
    [fr, dfr, Qm, Rm, C0, Rs] = c2v(coeffs);
    
    Z = model(n, exp(logf), fr, dfr, Qm, Rm, C0, Rs);
    logZ = log(abs(Z));
    %logZ = abs(Z);
    
    if debug_plot
        fig bvdfit:coeffs;
        lineappend(coeffs);

        fig bvdfit:fitdebug;
        plot(logf, logZ, 'LineWidth', 0.5);
        recolor winter;
        drawnow;
    end    
end

function Z = model(n, f, fr, dfr, Qm, Rm, C0, Rs)
    N  = round((n-1)/2);
    fr = (fr-N*dfr):dfr:(fr+N*dfr);
    
    n = length(fr);
    m = length(f);

    Cm = 1./(2*pi*fr.*Qm.*Rm);
    Lm = (Qm.*Rm)./(2*pi*fr);
    Rm = repmat(Rm, m, n);
    
    ZCm = 1./(2i*pi*f*Cm);
    ZLm =     2i*pi*f*Lm;
    ZC0 = 1./(2i*pi*f*C0);
    
    Zm = 1./sum(1./(ZLm + ZCm + Rm), 2);
    Z0 = ZC0 + Rs;
    Z  = Zm.*Z0./(Zm + Z0);
    
%     fprintf('\nfitmodel:   \n')
%     fprintf('fr  = %.3fGHz \n', fr(N)/1e9);
%     fprintf('dfr = %.2f MHz\n', dfr/1e6);
%     fprintf('Qm  = %.0f    \n', Qm);
%     fprintf('Rm  = %.2f    \n', Rm(1));
%     fprintf('C0  = %.2f pF \n', C0*1e12);
%     fprintf('Rs  = %.2f    \n', Rs);
end

function initial = coeff_estimate(f, z, n)

    f = f(:);
    z = z(:);

    [fr, fa, ir, ia] = resonance_find(f, z, n);
    
    za = abs(z(ia)); a = (2*pi*fa).^2;
    zr = abs(z(ir)); r = (2*pi*fr).^2; 
    
    c0 = 1./(2*pi*f.*abs(z));
    C0 = median(c0, 'omitnan');
    
    Cm = (a-r)./r*C0;
    Lm = 1./(a-r)/C0;
    
    Qm = a./(a-r) ./ (2*pi*fr*C0.*zr);
    Rm = 1./Qm .* sqrt(Lm./Cm);
    Rs = 1;
    
    % pick the fr in the center of the frequency-range
    dfr = median(abs(diff(fr)));
    [~, i] = min(abs(median(fr) - fr));
    
    initial = v2c(fr(i), dfr, median(Qm), median(Rm), C0, Rs);
    
%     fig bvdfit:initial; cla; yscale log; 
%     plot(f, c0);
%     plot(f, sort(c0));
%     plot([min(f) max(f)], ones(1,2)*C0, 'k--');
end

function data = fill_data_struct(initial, final, n, f, zfit, zraw)
    [frf, dfrf, Qmf, Rmf, C0f, Rsf] = c2v(final);
    [fri, dfri, Qmi, Rmi, C0i, Rsi] = c2v(initial);
    
    Cm = 1./(2*pi*frf.*Qmf.*Rmf);
    Lm = 1./(2*pi*frf).*Qmf.*Rmf;
    
    model = @(f)exp(wrapper(n, final,log([f(:)])));
    
    fa = sqrt(1./(Lm.*Cm).*(1 + Cm./C0f))/(2*pi);
    %kt2 = Cm./(C0f+Cm);
    %kt2 = pi/2*sqrt(C0f/(Cm + C0f))/tan(pi/2*sqrt(C0f/(Cm + C0f)));
    kt2 = pi/2 * frf/fa / tan(pi/2 * frf/fa);% * C0f/Cm / 7.5;
    
    i = (f > frf-dfrf*.6 & f < frf+dfrf*.6);
    fi = f(i);
    [zmax, imax] = max(zfit(i));
    [zmin, imin] = min(zfit(i));
    fmax = fi(imax);
    fmin = fi(imin);
    
    data.initial.fr = fri;
    data.initial.dfr = dfri;
    data.initial.R = Rmi;
    data.initial.Q = Qmi;
    data.initial.C0 = C0i;
    data.initial.Rs = Rsi;
    data.fr   = frf;
    data.fa   = fa;
    data.zr   = model(frf);
    data.za   = model(fa);
    data.fmin = fmin;
    data.fmax = fmax;
    data.zmin = zmin;
    data.zmax = zmax;
    data.dfr  = dfrf;
    data.R    = Rmf;
    data.Q    = Qmf;
    data.Rm   = Rmf;
    data.Qm   = Qmf;
    data.C0   = C0f;
    data.Cm   = Cm;
    data.Lm   = Lm;
    data.Rs   = Rsf;
    data.kt2  = kt2;
    data.model = @(f)exp(wrapper(n, final,log([f(:)])));
    data.zfit = zfit;
    data.zraw = zraw;
    data.f   = f;
    data.coeff_final = final;
    data.valid = ~any(isnan(final)) ... 
        && all(frf > min(f)) && all(frf < max(f)) ... 
        && all([Qmf Rmf Cm Lm] > 0);% ...
       % && abs(frf-fri) < dfri/100;
end

function lineappend(x, y)
    ax = gca;
    lines = ax.Children;
    n = length(x);
    
    if nargin < 2
        y = x;
        x = [];
    end
    if length(y) ~= n
        error 'x and y values must be of equal length';
    end
    
    warning off MATLAB:gui:array:InvalidArrayShape;
    for k = 1:n
        if ~isempty(x)
            lines(k).XData = [lines(k).XData x(k)];
        else
            N = length(lines(k).XData);
            lines(k).XData = [lines(k).XData N+1];
        end
        lines(k).YData = [lines(k).YData y(k)];
    end
end

function prepareCoeffPlot()
    colors = get(gca, 'ColorOrder');
    
    line('DisplayName', 'f_r', 'Color', colors(1,:));
    line('DisplayName', '\Delta{}f_r', 'Color', colors(1,:));
    line('DisplayName', 'R_m', 'Color', colors(2,:));
    line('DisplayName', 'Q_m', 'Color', colors(3,:));
    line('DisplayName', 'C_0', 'Color', colors(4,:));
    line('DisplayName', 'R_s', 'Color', colors(4,:));
    legend show; legend location EastOutside;
    xscale lin;

end
