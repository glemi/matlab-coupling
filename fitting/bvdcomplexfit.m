function [data, zfit] = bvdcomplexfit(f, z, varargin)
    hf = gcf;
    global debug_plot;
    
    if mod(length(varargin),2) == 1
        weights = varargin{1};
        varargin(1) = [];
    else
        weights = ones(size(f));
    end
    weights = [weights(:); weights(:)];
    
    opts = statset(varargin{:});
    % initialize variables
    logf    = log(f);
    initial = coeff_estimate(f, abs(z)); 
    final   = [NaN NaN NaN NaN NaN];   

    if ~any(isnan(initial))
        if debug_plot
            fig bvdfit:fitdebug; clf;
            plot(logf, real(z), '--', 'LineWidth', 2);
            plot(logf, imag(z), '--', 'LineWidth', 2);

            fig bvdfit:coeffs; clf;
            colors = get(gca, 'ColorOrder');
            line('DisplayName', 'f_r', 'Color', colors(1,:));
            line('DisplayName', 'R_m', 'Color', colors(2,:));
            line('DisplayName', 'Q_m', 'Color', colors(3,:));
            line('DisplayName', 'C_0', 'Color', colors(4,:));
            line('DisplayName', 'R_s', 'Color', colors(5,:));
            legend show; legend location EastOutside;
            xscale lin;

            final = nlinfit([logf; logf], [real(z); imag(z)], @wrapper, initial, opts, 'Weights', weights);
            zfit = wrapper(final, [logf; logf]);
            n = length(zfit);
            zfit = complex(zfit(1:n/2), zfit(1+n/2:end));
            
            fig bvdfit:fitdebug;
            plot(logf, real(zfit), '-', 'LineWidth', 2);
            plot(logf, imag(zfit), '-', 'LineWidth', 2); 
        else
            final = nlinfit([logf; logf], [real(z); imag(z)], @wrapper, initial, opts, 'Weights', weights);
            zfit = wrapper(final, [logf; logf]);
            n = length(zfit);
            zfit = complex(zfit(1:n/2), zfit(1+n/2:end));
        end
        
        %zfit = exp(zfit);
        data = fill_data_struct(initial, final, f, zfit, z);
    else
        zfit = nan(size(f));
        data = fill_data_struct(initial, final, f, zfit, z);
    end
    
    
    figure(hf);
end

function coeffs = v2c(fr, R, Q, C0, Rs)
    coeffs = [fr R Q C0 Rs].*[1e-9 1 1e-3 1e12 1];
end

function [fr, R, Q, C0, Rs] = c2v(coeffs)
    fr = coeffs(1)*1e9;
    R  = coeffs(2); 
    Q  = coeffs(3)*1e3;
    C0 = coeffs(4)*1e-12;
    Rs = abs(coeffs(5));
end

function Zri = wrapper(coeffs, logf)
    global debug_plot;
    
    coeffs = real(coeffs);
    fr = coeffs(1)*1e9;
    R  = coeffs(2);
    Q  = coeffs(3)*1e3;
    C0 = coeffs(4)*1e-12;
    Rs = coeffs(5);
    
    n = length(logf);
    logf = logf(1:n/2);
    Z = model(exp(logf), fr, R, Q, C0, Rs);
    
    Zri = [real(Z); imag(Z)];
     
    if debug_plot
        fig bvdfit:coeffs;
        lineappend(coeffs);

        fig bvdfit:fitdebug;
        plot(logf, real(Z), 'LineWidth', 0.5);
        plot(logf, imag(Z), 'LineWidth', 0.5);
        recolor winter;
        drawnow;
    end    
end

function Z = wrapper1(coeffs, f)
    coeffs = real(coeffs);
    
    fr = coeffs(1)*1e9;
    R  = coeffs(2);
    Q  = coeffs(3)*1e3;
    C0 = coeffs(4)*1e-12;
    Rs = coeffs(5);
    
    Z = model(f, fr, R, Q, C0, Rs);
end

function Z = model(f, fr, Rm, Qm, C0, R0)

    Cm = 1/(2*pi*fr*Qm*Rm);
    Lm = (Qm*Rm)/(2*pi*fr);
    
    ZCm = 1./(2i*pi*f*Cm);
    ZLm =     2i*pi*f*Lm;
    ZC0 = 1./(2i*pi*f*C0);
    
    Zm = ZLm + ZCm + Rm;
    Z0 = ZC0 + R0;
    Z  = Zm.*Z0./(Zm + Z0);% + abs(Rs);
end


function initial = coeff_estimate(f, z)

    [fr, fa, ir, ia] = resonance_find(f, z);
    za = abs(z(ia(1))); a = (2*pi*fa)^2;
    zr = abs(z(ir(1))); r = (2*pi*fr)^2; 
    
    c0 = 1./(2*pi*f.*abs(z));
    C0 = median(c0, 'omitnan');
    
    Cm = (a-r)/r*C0;
    Lm = 1/(a-r)/C0;
    
    Qm = a/(a-r) / (2*pi*fr*C0*zr);
    Rm = 1/Qm * sqrt(Lm/Cm);
    Rs = 5;
    
    initial = v2c(fr, Rm, Qm, C0, Rs);
    
%     fig bvdfit:initial; cla; yscale log; 
%     plot(f, c0);
%     plot(f, sort(c0));
%     plot([min(f) max(f)], ones(1,2)*C0, 'k--');
end

function data = fill_data_struct(initial, final, f, zfit, zraw)
    [frf, Rmf, Qmf, C0f, Rsf] = c2v(final);
    [fri, Rmi, Qmi, C0i, Rsi] = c2v(initial);
    
    Cm = 1/(2*pi*frf*Qmf*Rmf);
    Lm = 1/(2*pi*frf)*Qmf*Rmf;
    
    fa = sqrt(1/(Lm*Cm)*(1 + Cm/C0f))/(2*pi);

    %kt2 = Cm/(C0f+Cm);
    kt2 = pi/2*sqrt(C0f/(Cm + C0f))/tan(pi/2*sqrt(C0f/(Cm + C0f)));
    kt2 = pi/2 * frf/fa / tan(pi/2 * frf/fa);% * C0f/Cm / 7.5;
    
    data.initial.fr = fri;
    data.initial.R = Rmi;
    data.initial.Q = Qmi;
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
    data.model = @(f)exp(wrapper(final,log(f)));
    data.model1 = @(f)wrapper1(final,f);
    data.zfit = zfit;
    data.zraw = zraw;
    data.f   = f;
    data.coeff_final = final;
    data.valid = ~any(isnan(final)) ... 
        && frf > min(f) && frf < max(f) ... 
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
