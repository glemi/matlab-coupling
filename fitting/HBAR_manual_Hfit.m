function HBAR_manual_Hfit( )
    global debug_plot; 

    fig hbarui:test; clf;
    [f, Z] = loadMeasData;
    
   % mplot(f, Z, true);
    [Z C0] = HBAR_removeparasitics_v2(f, Z);
    %debug_plot = true;
    H = (2i*pi*f*C0).*Z;
    mplot(f, H, true);
    %debug_plot = false;

    hbarui = HBAR_paramui;
    hbarui.Callback = @(c)callback(hbarui, c);
end

function mplot(f, Z, keep)
    opt keep logical false;
    persistent hlines;
    
    curves = {'env:abs:upper' 'env:abs:lower' 'env:ph:upper' 'env:ph:lower' 'ripples:keff' 'ripples:Cm' 'ripples:diff(fr)' 'renv:abs:u/l' 'ripples:C0' 'ripples:Qm'};
    %curves = {'direct:abs' 'direct:ph'};
    
    
    [M, F] = HBAR_postprocess(f, Z, curves, 100, keep);
    fig hbarui:test;
    try delete(hlines); end;
    if keep
        h = HBAR_Mplot(F, M, curves, '-');
        delete(h(1:4)); 
        ax = subplot(4,2,1); ax.ColorOrderIndex = ax.ColorOrderIndex -2;
        %ax = subplot(2,1,1);
        plot(f/1e9, abs(Z));
        
        ax = subplot(4,2,2); ax.ColorOrderIndex = ax.ColorOrderIndex -2;
       % ax = subplot(2,1,2);
        plot(f/1e9, phase(Z)*180/pi);
    else
        hlines = HBAR_Mplot(F, M, curves, 'k-');
    end
end

function callback(hbarui, config)
    params = hbarui.params;
    HBAR_print('plain', params, config);
    
    f = [1e9:.1e6:4e9]';
    [Z, hbar] = HBAR_v3(f, config);
    mplot(f, hbar.H);
    
    names = {'kt2', 'C0'};
    values = [hbar.kt2 hbar.C0];
    latex = HBAR_print('blocklatex', names, values);
    latexfigure('hbarvalues', latex);
end


function [f, Z] = loadMeasData
    file = 'data\CTI_01_02_00C.s1p';
    data = read_s1p(file);
    f = data.f;
    Z = squeeze(data.z(1,1,:));
end