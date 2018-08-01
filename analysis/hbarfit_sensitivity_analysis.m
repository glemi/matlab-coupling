function  hbarfit_sensitivity_analysis(vary, xscale, xunit)
    yscale = 0.01; yunit = '%';

    f = 1e9:0.1e6:5e9; m = length(f);
    config = HBAR_loadconfig('HBAR_config.txt');
    [Zref, hbar_ref] = HBAR_v3(f, config);
    
    fitparams = {'ePiezo' 'cPiezo' 'epsPiezo' 'tPiezo' 'tBotEl' 'tTopEl' 'tSubst'};
    units     = {'C/m2'   'Pa'     '1'        'm'      'm'      'm'      'm'};
    scale     = [ 1        1e9      1          1e-9     1e-9     1e-9    1e-6];
    
    i = strcmp(fitparams, vary);
    fitparams(i) = []; units(i) = []; scale(i) = [];
    N = length(fitparams);
    
    x0 = HBAR_getconfigvalues(config, vary);
    y0 = HBAR_getconfigvalues(config, fitparams);
    kt20 = hbar_ref.kt2;
    
    xpercent = [80:5:120];
    values = xpercent*x0/100;
    
    i = strcmp(fitparams, vary);
    fitparams(i) = [];
    units(i) = [];

    n = length(values);
    
    Zfit = nan(m,n);
    y = nan(n, N);
    kt2 = nan(1, n);

    for k = 1:n
        try
            fixedparam  = {vary values(k)};
            [Zfit(:,k), fit] = hbarmultifit(f, Zref, config, fitparams, fixedparam, 'MaxIter', 5);
            fitconfig   = fit.hbar.config;
            y(k,:)      = HBAR_getconfigvalues(fitconfig, fitparams);
            kt2(k)      = fit.hbar.kt2;
        catch err
            message = 'error occurred while testing sensitivity to %s (iteration #%d)';
            warning(message, vary, k);
            errdisp(err);
        end
    end
    
    x = values;
    %y = [hbar.kt2];
    
    X = x/xscale;  
    X0 = x0/xscale;
    
    fig(sprintf('hbarsens:%s', vary)); clf;
    
    
    subplot(1,2,1);
        aleft = dualax('left'); 
        aleft.XLim = [X(1) X(end)];
        xlabel(sprintf('%s [%s]', vary, xunit));
        ylabel('k_t^2 [%]');

        plot([X(1) X(end)], [kt20 kt20]*100, '--');
        plot(X, kt2*100);
        
        axis manual;
        crosshair(X0, kt20*100);
        
        aright = dualax('right');
        aright.XLim = [xpercent(1) xpercent(end)];
        aright.XTickMode = 'auto';
        aright.YLim = aleft.YLim/kt20;
        aright.YTickMode = 'auto';
        

        xlabel '% of Reference Value';
        ylabel '% of Reference Value';
        htitle = title(sprintf('k_t^2 sensitivity on %s', vary));   
        
    subplot(1,2,2);
        n = size(y, 2);
        aleft = dualax('left');
        xlabel(sprintf('%s [%s]', vary, xunit));
        %ylabel '% of Reference Value';
        aleft.XLim = [X(1) X(end)];
        plot([X(1) X(end)], [100 100], 'k--');
        skiplegend; skipcolor;
        
        for k = 1:n
            ypercent = y(:,k)/y0(k)*100;
            plot(X, ypercent, 'DisplayName', fitparams{k});
        end
        legend show; legend location best;
        
        aright = dualax('right');
        crosshair(X0, 100);
        aright.XLim = [xpercent(1) xpercent(end)];
        aright.XTickMode = 'auto';
        aright.YTick = [];
        xlabel '% of Reference Value';
        axes(aleft);
end

