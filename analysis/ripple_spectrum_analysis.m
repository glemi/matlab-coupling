function ripple_spectrum_analysis(recompute)
    opt recompute boolean false;
    persistent Zripple Zclean f;

    if isempty(Zripple) || recompute
        nostack = {'tSi', 1e-9, 'tSiO', 1e-9, 'tBotE', 1e-9, 'tTopE', 1e-9};
        [Zripple, f] = rstor_model('Qpiezo', 10, 'fstep', 2e3);
        Zclean = rstor_model('fstep', 10e3, nostack{:});
    end
    
    z = abs(Zripple);
    
    fig spectrum:zel; clf;
    plot(f, z); xscale log; yscale log;
    
    [~, ia] = findpeaks((z));
    [~, ir] = findpeaks(1./(z));

    fr = f(ir);
    fa = f(ia);
    
    %[fr, fa] = pair_resonances(fr, fa);
    
    k2eff = (fa.^2 - fr.^2)./(fa.^2);
    k2eff_Cheeke = pi^2/4 * fr./fa .* (1 - fr./fa);
    k2eff_Naik = pi/2*fr./fa./tan(pi/2*fr./fa);
    
    fig spectrum:k2eff; clf;
    %plot(fr/1e9, k2eff*1000,        'DisplayName', 'Standard');
    %plot(fr/1e9, k2eff_Cheeke*1000, 'DisplayName', 'Cheeke');
    plot(fr/1e9, k2eff_Naik*1000,   'DisplayName', 'Naik');
    xlabel 'Frequency [GHz]';
    ylabel '[10^{-3}]';
    label([0.4 0.7], '$$k_{\rm eff}^2 = \frac{\pi}{2}\frac{f_r}{f_a}\tan{}\left( \frac{\pi}{2} \frac{f_r}{f_a} \right)$$');
    %label([0.35 0.15], '$$k_{\rm eff}^2 = \frac{f_a^2 - f_r^2}{f_a^2}$$');
    title 'Effective Coupling Factor'
    %legend show;
    
    fig spectrum:resdiff; clf;
    offset = 5.8e6;
    dfr = diff(fr);
    dfa = diff(fa);
    hlr = plot(fr(2:end)/1e9, (dfr -offset)/1e3);
    hla = plot(fa(2:end)/1e9, (dfa -offset)/1e3);
    hlr.DisplayName = 'Resonances: $\Delta{}f_r$';
    hla.DisplayName = 'Antiresonances: $\Delta{}f_a$';
    ylim([-100 100]);
    xlabel 'Frequency [GHz]';
    ylabel 'Spacing [kHz]';
    title 'Resonance/Antiresonance Spacing' 
    legend show; legend location south;
    set(legend, 'Interpreter', 'latex');
end

function [fr1, fa1] = pair_resonances(fr, fa)
    
    fr = fr/1e6;
    fa = fa/1e6;

    dfr = diff(fr);
    dfa = diff(fa);
    
    mdfr = median(dfr);
    mdfa = median(dfa);
    
    fig spectrum:resdiff; clf;
    
    plot(fr(2:end), dfr);
    plot(fa(2:end), dfa);
%     
    
%     mr = mod((fr-fr(1))/mdfr, 1);
%     ma = mod((fa-fa(1))/mdfa, 1); 
%     
%     plot(fr(2:end), diff(mr));
%     plot(fa(2:end), diff(ma));
    
%     plot(fr, ((mr)));
%     plot(fa, ((ma)));
    

    % Fr: columns, Fa: rows
    [Fr, Fa] = meshgrid(fr, fa);
    
    dFra = abs(Fr-Fa);
    %bdfra = (dFra < mdfr*0.4999999999999);
    bdfra = (dFra < mdfr*.5);
    
    fig spectrum:matrix; cla;
    spy(bdfra);
    axis tight;
    
    [i, j] = find(bdfra);
    
    
    fig spectrum:indices; clf;
    plot(diff(i));
    plot(diff(j));
    
    fr1 = fr(i)*1e6;
    fa1 = fa(j)*1e6;
end




