% function Zclean = HBAR_removeparasitics(f, Z)
function [Zclean, fit] = HBAR_removeparasitics_v4(f, Z, variables)
    
    zr = deripple(f, real(Z));
    zi = deripple(f, imag(Z));
    Zc = complex(zr, zi);
    
    % optional variables: R0 Rp Rc Lc U    
    if nargin < 3
        variables = {'Rc' 'Lc' 'R0' 'td'}; 
    end
    fit = bvdextfit_v2(f, Zc, variables);
    
    Zf = fit.Zfit;    
    
    fit.Lc = fit.Lc;
    Rc = fit.Rc;
    Rp = fit.Rp;
    R0 = fit.R0;
    
    ZC0 = 1./(2i*pi*f*fit.C0);
    ZCm = 1./(2i*pi*f*fit.Cm);
    ZLm =     2i*pi*f*fit.Lm;
    ZLc =     2i*pi*f*fit.Lc;
    Zx1  =    (fit.X1*2*pi*f);
    Zx2  =    (fit.X2*2*pi*f).^2;
    Zy1  = 1./(fit.Y1*2*pi*f);
    Zy2  = 1./(fit.Y1*2*pi*f).^2;
    Zm  = ZCm + ZLm + fit.Rm;
    
    function z1 = deembed(z, variable)
        switch variable
            case 'Rc',  z1 = z - Rc;
            case 'Lc',  z1 = z - ZLc;
            case 'Cs',  z1 = z - ZCs;
            case 'Rp',  z1 = Rp*z./(Rp - z);
            case 'X1',   z1 = z - Zx1;
            case 'X2',   z1 = z - Zx2;
            case 'Y1',   z1 = z - Zy1;
            case 'Y2',   z1 = z - Zy2;
            case 'R0'
                z1 = 1./(1./z - 1./Zm) - ZC0 - R0;
                z1 = 1./(1./(z1 + ZC0) + 1./Zm);
            case 'td',  z1 = z; % no way to extract this;
        end
    end

    n = length(variables);
    for k = 1:n
        var = variables{k};
        Z  = deembed(Z,  var);
        Zc = deembed(Zc, var);
        Zf = deembed(Zf, var);
    end
    
    Zclean = Z;
end
