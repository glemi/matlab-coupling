function [data, zfit] = bvdfit(f, z, n, varargin)
    hf = gcf;
    global debug_plot;
    
    opts = statset(varargin{:});
    % initialize variables
    logf    = log(f);
    initial = coeff_estimate(f, z, n); 
    final   = nan(size(initial));   

    if ~any(isnan(initial))
        if debug_plot
            fig bvdfit:fitdebug; clf;
            plot(logf, log(abs(z)));

            fig bvdfit:coeffs; clf;
            prepareCoeffPlot(n);

            final = nlinfit(logf, log(abs(z)), @wrapper, initial, opts);
            logz = wrapper(final, logf);

            fig bvdfit:fitdebug;
            plot(logf, abs(logz), '--', 'LineWidth', 2);        
        else
            final = nlinfit(logf, log(abs(z)), @wrapper, initial, opts);
            logz = wrapper(final, logf);
        end
        
        zfit = exp(logz);
        data = fill_data_struct(initial, final, f, zfit, z);
    else
        zfit = nan(size(f));
        data = fill_data_struct(initial, final, f, zfit, z);
    end
    
    
    figure(hf);
end

function coeffs = v2c(fr, Qm, Rm, C0, Rs)
    fr = fr*1e-9;
    Qm = Qm*1e-3;
    C0 = C0*1e12;
    coeffs = [fr(:)' Qm(:)' Rm(:)' C0 Rs];
end

function [fr, Qm, Rm, C0, Rs] = c2v(coeffs)
    n = (length(coeffs) - 2)/3;

    fr  = coeffs(0*n+1:1*n)*1e9;
    Qm  = coeffs(1*n+1:2*n)*1e3;
    Rm  = coeffs(2*n+1:3*n); 
    
    C0  = coeffs(3*n+1)*1e-12;
    Rs  = abs(coeffs(3*n+1));
end

function logZ = wrapper(coeffs, logf)
    global debug_plot;
    
    [fr, Qm, Rm, C0, Rs] = c2v(coeffs);
    
    Z = model(exp(logf), fr, Qm, Rm, C0, Rs);
    logZ = log(abs(Z));
    
    if debug_plot
%         fig bvdfit:coeffs;
%         lineappend(coeffs);

        fig bvdfit:fitdebug;
        plot(logf, logZ, 'LineWidth', 0.5);
        recolor winter;
        drawnow;
    end    
end

function Z = model(f, fr, Qm, Rm, C0, Rs)
    n = length(fr);
    m = length(f);
    
%     n = 1;
%     fr = fr(1);
%     Rm = Rm(1);
%     Qm = Qm(1);

    Cm = 1./(2*pi*fr.*Qm.*Rm);
    Lm = (Qm.*Rm)./(2*pi*fr);
    Rm = repmat(Rm, m, 1);
    
    ZCm = 1./(2i*pi*f*Cm);
    ZLm =     2i*pi*f*Lm;
    ZC0 = 1./(2i*pi*f*C0);
    
    Zm = 1./sum(1./(ZLm + ZCm + Rm), 2);
    Z0 = ZC0; %+ R0;
    Z  = Zm.*Z0./(Zm + Z0) + abs(Rs);
end

function initial = coeff_estimate(f, z, n)

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
    
    initial = v2c(fr, Qm, Rm, C0, Rs);
    
%     fig bvdfit:initial; cla; yscale log; 
%     plot(f, c0);
%     plot(f, sort(c0));
%     plot([min(f) max(f)], ones(1,2)*C0, 'k--');
end

function data = fill_data_struct(initial, final, f, zfit, zraw)
    [frf, Rmf, Qmf, C0f, Rsf] = c2v(final);
    [fri, Rmi, Qmi, C0i, Rsi] = c2v(initial);
    
    n = length(frf);
    i = round((n+1)/2);
    frf = frf(i);
    Rmf = Rmf(i);
    Qmf = Qmf(i);
    
    
    Cm = 1./(2*pi*frf.*Qmf.*Rmf);
    Lm = 1./(2*pi*frf).*Qmf.*Rmf;
    
    fa = sqrt(1./(Lm.*Cm).*(1 + Cm./C0f))/(2*pi);
    kt2 = Cm./(C0f+Cm);
    
    data.initial.fr = fri(i);
    data.initial.R = Rmi(i);
    data.initial.Q = Qmi(i);
    data.initial.C0 = C0i;
    data.initial.Rs = Rsi;
    data.fr  = frf;
    data.fa  = fa;
    data.R   = Rmf;
    data.Q   = Qmf;
    data.C0  = C0f;
    data.Cm  = Cm;
    data.Lm  = Lm;
    data.Rs  = Rsf;
    data.kt2 = kt2;
    data.model = @(f)exp(wrapper(final,log([f(:)])));
    data.zfit = zfit;
    data.zraw = zraw;
    data.f   = f;
    data.coeff_final = final;
    data.valid = ~any(isnan(final)) ... 
        && all(frf > min(f)) && all(frf < max(f)) ... 
        && all([Qmf Rmf Cm Lm] > 0);
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

function prepareCoeffPlot(n)
    colors = get(gca, 'ColorOrder');
    
    for k = 1:n 
        line('DisplayName', sprintf('f_r %d', k), 'Color', colors(1,:));
        line('DisplayName', sprintf('R_m %d', k), 'Color', colors(2,:));
        line('DisplayName', sprintf('Q_m %d', k), 'Color', colors(3,:));
    end
    line('DisplayName', 'C_0', 'Color', colors(4,:));
    line('DisplayName', 'R_s', 'Color', colors(4,:));
    legend show; legend location EastOutside;
    xscale lin;

end
