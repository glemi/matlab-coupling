function ripple_analysis(parameter, values, parName, parUnit, parScale, axScale)
    opt parUnit char '';
    opt parScale double 1;
    opt axScale char lin;
    opt parName char parameter;

    f = (1e9:.1e6:5e9)';
    n = length(values);
    values = values*parScale;
    
    % Default values
    tPiezo = 800e-9;
    tSubst = 725e-6; 
    tOxide = 200e-9;
    tBotEl = 200e-9;
    tTopEl = 100e-9;
    
    for k = 1:n 
        arg_var = {parameter, values(k)};
        
        arg_PS  = {'tPiezo', tPiezo, 'tSi', tSubst, 'tSiO', tOxide, 'tBotE', 0,      'tTopE', 0  };
        arg_PSE = {'tPiezo', tPiezo, 'tSi', tSubst, 'tSiO', tOxide, 'tBotE', tBotEl, 'tTopE', tTopEl};

        [Z_PS,  ~,  par_FPS]  = HBAR(f, arg_PS{:},  arg_var{:});
        [Z_PSE, ~, par_FPSE]  = HBAR(f, arg_PSE{:}, arg_var{:});
        
        min_dfr = par_FPSE.frSi*0.9;
%         ps_ripples = multiRippleFit(f, Z_PS, 3, 5e6);
%         ps_data(k).ripples = ps_ripples; 
%         ps_data(k).param = par_FPS;
%         ps_data(k).zel = Z_PS;

        pse_ripples = multiRipplefit(f, abs(Z_PSE), 3, min_dfr);
        pse_data(k).ripples = pse_ripples;
        pse_data(k).param = par_FPSE;
        pse_data(k).zel = Z_PSE;
        
        if ~isempty(parUnit)
            names{k} = sprintf('%.1f %s', values(k)/parScale, parUnit);
        else
            names{k} = sprintf('%s', scinot(values(k)/parScale));
        end
    end
    
    if strcmp(parName, 'kt')
        par = [pse_data.param];
        values = [par.k2];
        parName = 'k_t^2';
        parScale = 0.01;
        parUnit = '%';
    end
    
    fname = sprintf('ripples:%s', parameter);
    fig(fname); clf; suptitle 'Extracted k_t^2';
    
    subplot(3,2,1); title '\Delta f_r'
    %dfrplot(f, ps_data,  values, parName, parUnit, parScale, axScale); 
    dfrplot(pse_data, values, parName, parUnit, parScale, axScale); 
    
    subplot(3,2,2); title 'k_{eff}^2 [%]'
    %keffplot(f, z, ps_data, values, parName, parUnit, parScale, axScale); 
    keffplot(pse_data, values, parName, parUnit, parScale, axScale); 
    
    %subplot(2,2,3); title 'k_t^2'
    %keffplot(f, z, ps_data, values, parName, parUnit, parScale, axScale); 
    ktplot(pse_data, values, parName, parUnit, parScale, axScale); 
    
    suptitle(['kt Extract vs ' parName]);
end

function dfrplot(data, values, parName, parUnit, parScale, axScale)    
    n = length(data);
    for k = 1:n
        label = sprintf('%.1f %s', values(k)/parScale, parUnit);
        r = data(k).ripples;
        
        fr  = [r.fr];
        dfr = [r.dfr];
        
        plot(fr/1e9, dfr/1e6, 'DisplayName', label);
        dfrdata = dfr_analyze(fr, dfr);
        dissolve(dfrdata); skipcolor;
        plot([f0min f0max]/1e9, [dfrmin dfrmax]/1e6, 'o');
        skiplegend; fillmarkers; 
    end
    legend show; legend location best;
end

function keffplot(data, values, parName, parUnit, parScale, axScale)    
    n = length(data);
    for k = 1:n
        label = sprintf('%.1f %s', values(k)/parScale, parUnit);
        r = data(k).ripples;
        
        fr  = [r.fr];
        keff = [r.kt2];
        plot(fr/1e9, keff*100, 'DisplayName', label);   
    end
    legend show; legend location best;
end

function ktplot(data, values, parName, parUnit, parScale, axScale)    
    n = length(data);
    for k = 1:n
        r = data(k).ripples;
        fr   = [r.fr];
        dfr  = [r.dfr];
        keff = [r.kt2];
        
        dfrdata(k)  = dfr_analyze(fr, dfr);
        cheeke(k)   = cheeke_compute(data(k).param);
        
        keff0       = interp1(fr, keff, dfrdata(k).f0min);
        kt2(k)      = cheeke(k).kt2Nfactor*keff0;
        kt2ref(k)   = data(k).param.k2;
    end
    

    subplot(3,2,3); title 'min & max of \Delta{}f_r spacing [MHz]';
    plot(values/parScale, [dfrdata.dfrmin]/1e6, '-o', 'DisplayName', 'min\{\Delta{}f_r\}'); 
    plot(values/parScale, [dfrdata.dfrmax]/1e6, '-o', 'DisplayName', 'max\{\Delta{}f_r\}');
    fillmarkers;
    legend show; legend location best;
    xlabel(sprintf('%s [%s]', parName, parUnit));
    
    subplot(3,2,4); title 'location of \Delta{}f_r min & max [GHz]';
    plot(values/parScale, [dfrdata.f0min]/1e9, '-o', 'DisplayName', 'min\{\Delta{}f_r\}'); 
    plot(values/parScale, [dfrdata.f0max]/1e9, '-o', 'DisplayName', 'max\{\Delta{}f_r\}');
    fillmarkers;
    legend show; legend location best;
    xlabel(sprintf('%s [%s]', parName, parUnit));
    
    subplot(3,2,5); title 'extracted k_t^2 [%]';
    plot(values/parScale, kt2*100, '-o', 'DisplayName', 'Cheeke''s method'); 
    plot(values/parScale, kt2ref*100, 'k--', 'DisplayName', 'reference value');
    fillmarkers;
    legend show; legend location best;
    xlabel(sprintf('%s [%s]', parName, parUnit));
    
    %axis equal; axis square;
end

function subax(item)
    switch item
        case 'P',   subplot(2,2,1); title 'Piezolayer Only (P)';
        case 'PE',  subplot(2,2,2); title 'Piezolayer + Electrodes (PE)';
        case 'PS',  subplot(2,2,3); title 'Piezolayer + Substrate (PS)';
        case 'PSE', subplot(2,2,4); title 'Full Stack (PSE)';
    end
end

