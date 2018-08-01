function stackcompare
    fig stackcomp:eins; clf; suptitle 'Impedance Spectra';
    fig stackcomp:zwei; clf; suptitle 'Extracted k_t^2';

    param = 'E33'; parName = 'e33'; values = [1.5:0.1:2.5];
    
    n = length(values);
    for k = 1:n 
        execute(param, values(k), parName);
        drawnow;
    end
end


function execute(parameter, value, parName, parUnit, parScale, axScale)
    opt parScale double 1;
    opt axScale char lin;
    opt parName char parameter;

    f = 1e9:1e6:10e9;
    arg_var = {parameter, value};
    
    arg_P   = {'tPiezo', 800e-9, 'tSi', 0,      'tSiO', 0,      'tBotE', 0,      'tTopE', 0};
    arg_PE  = {'tPiezo', 800e-9, 'tSi', 0,      'tSiO', 0,      'tBotE', 200e-9, 'tTopE', 100e-9};
    arg_PS  = {'tPiezo', 800e-9, 'tSi', 725e-6, 'tSiO', 200e-9, 'tBotE', 0,      'tTopE', 0  };
    arg_PSE = {'tPiezo', 800e-9, 'tSi', 725e-6, 'tSiO', 200e-9, 'tBotE', 200e-9, 'tTopE', 100e-9};
    
    [Z_P,   Z_FP,   par_FP   ]  = HBAR(f, arg_P{:}, arg_var{:});
    [Z_PE,  Z_FPE,  par_FPE  ]  = HBAR(f, arg_PE{:}, arg_var{:});
    [Z_PS,  Z_FPS,  par_FPS  ]  = HBAR(f, arg_PS{:}, arg_var{:});
    [Z_PSE, Z_FPSE, par_FPSE ]  = HBAR(f, arg_PSE{:}, arg_var{:});
    
    fig stackcomp:eins; 
    subax P;   zplot(f, Z_P);    zplot(f, Z_FP);  
    subax PE;  zplot(f, Z_PE);   zplot(f, Z_FPE);
    subax PS;  zplot(f, Z_PS);   zplot(f, Z_FPS);
    subax PSE; zplot(f, Z_PSE);  zplot(f, Z_FPSE);
    
    fig stackcomp:zwei; 
    subax P;   ktplot(f, Z_P,   par_FP  ); % ktplot(f, Z_FP,   par_FP  );  
    subax PE;  ktplot(f, Z_PE,  par_FPE ); % ktplot(f, Z_FPE,  par_FPE );
    subax PS;  ktplot(f, Z_PS,  par_FPS ); % ktplot(f, Z_FPS,  par_FPS );
    subax PSE; ktplot(f, Z_PSE, par_FPSE); % ktplot(f, Z_FPSE, par_FPSE);
end

function zplot(f, z)
    plot(f/1e9, abs(z));
    xscale log; yscale log;
    xlabel 'frequency [GHz]';
    ylabel 'Z_{el} [\Omega]';
end

function ktplot(f, z, par)
    kt2_ref = par.k2;
    data = bvdfit(f', abs(z));
    growplot(kt2_ref*100, data.kt2*100, '-o');
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
