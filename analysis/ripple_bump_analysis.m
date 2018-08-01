function ripple_bump_analysis( )

    [f, Z] = loadMeasData;
    
    fit = HBAR_ripplefit(f, abs(Z));
    
    fig ripple_bump:eins; clf;
    h1 = subplot(2,2,1);
    
    Zfit = NaN(size(f));
    
    n = length(fit);
    for k = 1:n
        x = fit(k);
        i = f > x.fr-x.dfr/2 & f < x.fr+x.dfr;
        Zfit(i) = x.model(f(i));
        
        fQm = getFixedQModel(x, min(1300, x.Q));
        ZfQ(i) = abs(fQm(f(i)));
    end
    
    plot(f, abs(Z));
    plot(f, Zfit);
    
    xscale log; 
    yscale log;
    
    h3 = subplot(2,2,3);
    Zd1 = abs(Z(:)) - Zfit(:);
    plot(f, Zd1);
    
    
    h2 = subplot(2,2,2);
    plot(f, abs(Z));
    plot(f, ZfQ(:));
    
    xscale log; 
    yscale log;
    
    h4 = subplot(2,2,4);
    Zd2 = abs(Z(:))-ZfQ(:);
    plot(f, Zd2);
    crosshair;
    
    linkaxes([h1 h3], 'x');
    linkaxes([h2 h4], 'x');
    
    fig ripple_bump:zwei;
    subplot(2,1,1);
    plot(f, Zd1);
    plot(f, Zd2);
    
    subplot(2,1,2);
    plot(f, Zd2 - Zd1)
    
    
    
%     zr = deripple(f, real(Z));
%     zi = deripple(f, imag(Z));
%     Zc = complex(zr, zi);
%     fit = bvdextfit(f, Zc);
%     Zf = fit.Zfit;    
%     
%     fig parasitics:null; clf;
%     rawdataplot(f, Z);
%     rawdataplot(f, Zc);
%     rawdataplot(f, Zf);
%         
%     Rp = fit.Rp;
%     Rc = fit.Rc;
%     
%     Z1 = Z-Rc;
%     Z2 = Rp*Z1./(Rp - Z1);
%     
%     fig parasitics:eins; clf;
%     rawdataplot(f, Z1);
%     
%     fig parasitics:zwei; clf;
%     rawdataplot(f, Z2);
    
end

function model = getFixedQModel(fit, Q)
    
    model = @(f)fixedQModel(51, f, fit.fr, fit.dfr, Q, fit.Rm, fit.C0, fit.Rs);
end


function Z = fixedQModel(n, f, fr, dfr, Qm, Rm, C0, Rs)
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
end

function [f, Z] = loadMeasData
    file = 'data\CTI_01_02_16C.s1p';
    data = read_s1p(file);
    f = data.f;
    Z = squeeze(data.z(1,1,:));
end