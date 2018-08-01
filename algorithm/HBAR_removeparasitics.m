
% function Zclean = HBAR_removeparasitics(f, Z)
function Zclean = HBAR_removeparasitics(f, Z)
    
    zr = deripple(f, real(Z));
    zi = deripple(f, imag(Z));
    Zc = complex(zr, zi);
    
    fit = bvdextfit(f, Zc);
    weights = setweights(f, fit.Zfit);
    fit = bvdextfit(f, Zc, 'Weights', weights);
    
    fr = fit.fr;
    Rm = fit.Rm;
    Qm = fit.Qm;
    C0 = fit.C0;
    
    Zideal = bvdidealmodel(f, fr, Rm, Qm, C0);
    Znoise = Zc - Zideal;
    Zclean = Z - Znoise;
end

function w = setweights(f, z)
    [fr, fa] = resonance_find(f, abs(z));
    df = fa-fr;
    
    i = f > fr-df/3 & f < fa+df;
    
    n = length(z);
    w = ones(n, 1);
    w(i) = n;
end

function Z = bvdidealmodel(f, fr, Rm, Qm, C0)
    Cm = 1/(2*pi*fr*Qm*Rm);
    Lm = (Qm*Rm)/(2*pi*fr);
    
    ZC0 = 1./(2i*pi*f*C0);
    ZCm = 1./(2i*pi*f*Cm);
    ZLm =     2i*pi*f*Lm;
    
    Zm = ZLm + ZCm + Rm;
    Z0 = ZC0;
    
    Z = Zm.*Z0./(Zm + Z0);
end
