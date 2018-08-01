function kt_extract_analysis(parameter, values, parName, parUnit, parScale, axScale)
    opt parUnit char '';
    opt parScale double 1;
    opt axScale char lin;
    opt parName char parameter;

    f = 1e9:.1e6:20e9;
    n = length(values);
    values = values*parScale;
    
    % Default values
    tPiezo = 800e-9;
    tSubst = 725e-6; 
    tOxide = 200e-9;
    tBotEl = 300e-9;
    tTopEl = 100e-9;
    
    for k = 1:n 
        arg_var = {parameter, values(k)};
        
        arg_P   = {'tPiezo', tPiezo, 'tSi', 0,      'tSiO', 0,      'tBotE', 0,      'tTopE', 0};
        arg_PE  = {'tPiezo', tPiezo, 'tSi', 0,      'tSiO', 0,      'tBotE', tBotEl, 'tTopE', tTopEl};
        arg_PS  = {'tPiezo', tPiezo, 'tSi', tSubst, 'tSiO', tOxide, 'tBotE', 0,      'tTopE', 0  };
        arg_PSE = {'tPiezo', tPiezo, 'tSi', tSubst, 'tSiO', tOxide, 'tBotE', tBotEl, 'tTopE', tTopEl};

        [Z_P(:,k),   Z_FP(:,k),   par_FP(k)  ]  = HBAR(f, arg_P{:},   arg_var{:});
        [Z_PE(:,k),  Z_FPE(:,k),  par_FPE(k) ]  = HBAR(f, arg_PE{:},  arg_var{:});
        [Z_PS(:,k),  Z_FPS(:,k),  par_FPS(k) ]  = HBAR(f, arg_PS{:},  arg_var{:});
        [Z_PSE(:,k), Z_FPSE(:,k), par_FPSE(k)]  = HBAR(f, arg_PSE{:}, arg_var{:});
        
        if ~isempty(parUnit)
            names{k} = sprintf('%.1f %s', values(k)/parScale, parUnit);
        else
            names{k} = sprintf('%s', scinot(values(k)/parScale));
        end
    end

    fname = sprintf('ktextract:%s', parameter);
    fig(fname); clf; suptitle 'Extracted k_t^2';
    subax P;   ktplot(f, Z_P,   par_FP,  values, parName, parUnit, parScale, axScale); 
    subax PE;  ktplot(f, Z_PE,  par_FPE, values, parName, parUnit, parScale, axScale); 
    subax PS;  ktplot(f, Z_PS,  par_FPS, values, parName, parUnit, parScale, axScale); 
    subax PSE; ktplot(f, Z_PSE, par_FPSE,values, parName, parUnit, parScale, axScale); 
    
    suptitle(['kt Extract vs ' parName]);
end

function ktplot(f, z, data, values, parName, parUnit, parScale, axScale)    
    n = size(z, 2);
    
    kt2_ref = [data.k2];
    plot(values/parScale, kt2_ref*100, '--');
    plot(NaN, NaN, '-o'); fillmarkers;
    ylim([min(kt2_ref)*0.5 max(kt2_ref)*1.5]*100);
    xlabel(sprintf('%s [%s]', parName, parUnit));
    xscale(axScale);
    drawnow;
    
    for k = 1:n
        fit(k) = bvdfit(f', abs(z(:,k)));
        growplot([values(k)/parScale NaN], [fit(k).kt2*100 NaN]);
        drawnow;
    end
    %axis equal; axis square;
end

function ktplot1(f, z, par)    
    n = size(z, 2);
    for k = 1:n
        kt2_ref(k) = par(k).k2;
        data(k) = bvdfit(f', abs(z(:,k)));
    end
    plot(kt2_ref*100, [data.kt2]*100, '-o');
    plot(kt2_ref*100, kt2_ref*100, '--');
    
    fillmarkers;
    axis equal; axis square;
end

function subax(item)
    switch item
        case 'P',   subplot(2,2,1); title 'Piezolayer Only (P)';
        case 'PE',  subplot(2,2,2); title 'Piezolayer + Electrodes (PE)';
        case 'PS',  subplot(2,2,3); title 'Piezolayer + Substrate (PS)';
        case 'PSE', subplot(2,2,4); title 'Full Stack (PSE)';
    end
end
