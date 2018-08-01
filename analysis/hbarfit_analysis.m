function hbarfit_analysis
    
    f = (0.5e9:.1e6:5e9);
    config = 'HBAR_config.txt';
    
    
    %[Z, ~, hbardata] = HBAR(f, 'E33', 3, 'tPiezo', 850e-9, 'tBotE', 300e-9, 'tTopE', 200e-9);
    %[Z, hbardata] = HBAR_v2(f, 'E33', 3, 'tPiezo', 850e-9, 'tBotE', 300e-9, 'tTopE', 200e-9);
    [Z, hbardata] = HBAR_v3(f, config,  'ePiezo', 3, 'tPiezo', 850e-9, 'tBotEl', 300e-9, 'tTopEl', 200e-9);
    z = deripple(f, Z);
    
    fig hbarfit:eins; clf;
    subplot(2, 2, 1);
    xscale log; yscale log;
    
    plot(f, abs(Z));
    plot(f, z);
    xlim([f(1) f(end)]);
    %axis manual;
    
    %[zfit, fitdata] = hbarfit(f, z);
    %[zfit, fitdata] = hbarfit_v2a(f, Z);
    fitparams = {'ePiezo' 'tPiezo' 'tBotEl' 'tTopEl'};
    units = {'C/m2' 'm' 'm' 'm'};
    %[zfit, fitdata] = hbarfit_v3(f, Z, config, fitparams);
    [zfit, fitdata] = hbarmultifit(f, Z, config, fitparams);
    
    fig hbarfit:eins; 
    plot(f, abs(zfit), 'LineWidth', 3);
    
    xscale log;
    yscale log;
    
    fprintf('actual:    %.1f%%\n', hbardata.kt2*100);
    fprintf('extracted: %.1f%%\n', fitdata.hbar.kt2*100);
    
    ltext = printdata(fitdata, units);
    ltext = [ltext sprintf('%8s: %.1f%%', 'k_t^2', fitdata.hbar.kt2*100)];
    label([0.05 0.2], sprintf('%s\n', ltext{:}), 'Times12');
end

function text = printdata(fitdata, units)
    n = length(fitdata.final);
    
    for k = 1:n
        number  = siPrefix(fitdata.final(k), units{k});
        text{k} = sprintf('%8s: %s', fitdata.names{k}, number);
    end
end