% function Zclean = HBAR_removeparasitics(f, Z)
function [Zclean, fit] = HBAR_removeparasitics_v2(f, Z)
    
    zr = deripple(f, real(Z));
    zi = deripple(f, imag(Z));
    Zc = complex(zr, zi);
    fit = bvdextfit(f, Zc);
    Zf = fit.Zfit;    
    
        
    Rp = fit.Rp;
    Rc = fit.Rc;
    Lc = fit.Lc;
    R0 = fit.R0;
    
    Z1 = Z-Rc-2i*pi*f*Lc -1./(R0*2*pi*f).^2;
    Z2 = Rp*Z1./(Rp - Z1);
    
    Zclean = Z2;
    C0 = fit.C0;
    
%     fig parasitics:null; clf;
%     rawdataplot(f, Z);
%     rawdataplot(f, Zc);
%     rawdataplot(f, Zf);
%     subplot(2,2,1);
%     label([0.6 0.7], sprintf('$C_0 = %.2f$\\,pF', C0*1e12));
%     
%     fig parasitics:eins; clf;
%     rawdataplot(f, Z1);
%     rawdataplot(f, Zc);
%     subplot(2,2,1);
%     label([0.6 0.7], sprintf('$R_c = %.2f\\,\\Omega$', Rc));
    
%     fig parasitics:zwei; clf;
%     rawdataplot(f, Z2);
%     rawdataplot(f, Zc);
%     subplot(2,2,1);
%     label([0.5 0.8], sprintf('$R_p = %.0f\\,\\rm M\\Omega$', Rp/1e6));
    
end
