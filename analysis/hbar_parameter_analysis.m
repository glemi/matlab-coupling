function hbar_parameter_analysis(param)

    f = [1e9:0.1e6:5e9]'; m = length(f);
    config = HBAR_loadconfig('HBAR_config.txt');
    
    x0 = HBAR_parameters(config, param);
    xpercent = [80:10:120];
    x = xpercent*x0/100;
    n = length(x);
    
    codes = {'ripples:C0' 'ripples:Cm' 'ripples:Lm' 'renv:abs:u/l' ... 
        'ripples:keff' 'ripples:diff(fr)'};
    
    %fig ripplefit:coeffs; clf;
    figname = sprintf('hbarpar:%s', param);
    fig(figname); clf;
    for k = 1:n
        [Z, hbar] = HBAR_v4(f, config, param, x(k));
        [M, F] = HBAR_postprocess(f, Z, codes, 400);
        fig(figname);
        HBAR_Mplot(F, M, codes); drawnow;
    end
    
    entries = HBAR_print('latex', param, x);
    L = legend(entries);
    L.Interpreter = 'latex';
    
    
end
