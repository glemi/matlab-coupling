function model_output_analysis(parameter, values, parName, parUnit, parScale, axScale)
    opt parUnit char '';
    opt parScale double 1;
    opt axScale char lin;
    opt parName char parameter;

    f = 1e9:.1e6:10e9;
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
    
    fname = sprintf('modeloutput:%s', parameter);
    fig(fname); clf;
    subax P;   zplot(f, Z_P);   % zplot(f, Z_FP);  
    subax PE;  zplot(f, Z_PE);  % zplot(f, Z_FPE);
    subax PS;  zplot(f, Z_PS);  % zplot(f, Z_FPS);
    subax PSE; zplot(f, Z_PSE); % zplot(f, Z_FPSE);
    
    subax P;
    legend(names); legend location best;
    
    suptitle(['Impedance Spectra vs ' parName]);
end


function zplot(f, z)
    plot(f/1e9, abs(z));
    xscale log; yscale log;
    xlabel 'frequency [GHz]';
    ylabel 'Z_{el} [\Omega]';
end


function subax(item)
    switch item
        case 'P',   subplot(2,2,1); title 'Piezolayer Only (P)';
        case 'PE',  subplot(2,2,2); title 'Piezolayer + Electrodes (PE)';
        case 'PS',  subplot(2,2,3); title 'Piezolayer + Substrate (PS)';
        case 'PSE', subplot(2,2,4); title 'Full Stack (PSE)';
    end
end
