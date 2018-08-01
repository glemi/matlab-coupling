function kt_extract_compare_models(varargin)

    if nargin >= 1 && ischar(varargin{1}) && strcmpi(varargin{1}, 'fbar')
        data = fbar_compute();
        
        fig ktext:fbar; clf; ktextplot(data); title 'FBAR';
        fig frext:fbar; clf; frplot(data);    title 'FBAR';
        fig fit:fbar;  clf;  rawplot(data);   title 'FBAR';
        
    else
        data = hbar_compute(varargin{:});
        pstr = par2str(varargin);
        fig(['ktext:' pstr]); clf; ktextplot(data); title(pstr);
        fig(['frext:' pstr]); clf; frplot(data);    title(pstr);
        fig(['fit:' pstr]); clf; rawplot(data);   title(pstr);
    end
end

function rawplot(data)
    n = length(data);
    i = floor([n n/2 1]);
    for k = 1:3
        item = data(i(k));
        
        fraw = item.hbar.f;
        Zraw = item.hbar.Z;
        ffit = item.bvd.f;
        Zfit = item.bvd.zfit;
        Zctr = item.Zc;
        
        subplot(3,1,k);
        plot(fraw/1e9, abs(Zraw), 'Displayname', '|Z_{raw}|');
        plot(ffit/1e9, Zctr, 'Displayname', '|Z_{avg}|');
        plot(ffit/1e9, Zfit, '--', 'Displayname', '|Z_{fit}|' );
        
        fpeak = [item.fr_peak NaN item.fa_peak];
        Zpeak = [item.Zr_peak NaN item.Za_peak];
        fbvd  = [item.fr_bvd  NaN item.fa_bvd];
        Zbvd  = [item.Zr_bvd  NaN item.Za_bvd];
        
        set(gca, 'ColorOrderIndex', 2);
        plot(fpeak/1e9, Zpeak, 'o', 'DisplayName', 'f_a, f_r (peak)');
        plot(fbvd/1e9, Zbvd,   'o', 'DisplayName', 'f_a, f_r (BvD)');
        fillmarkers;
        
        label([0.75 0.8], sprintf('$k_t^2 = %.1f\\,\\%%$', item.kt2_ref*100));
        
        xscale log;
        yscale log;
        xlim(fraw([1 end])/1e9);
    end
    xlabel 'frequency [GHz]';
    
    subplot(3,1,1); axis manual;
    legend show; legend location SouthWest;
end

function frplot(data)
    x = real([data.kt2_ref]*100);
    
    plot(x, [data.fr_bvd]/1e9, 'o-', 'DisplayName', 'f_r (bvd)');
    plot(x, [data.fa_bvd]/1e9, 'o-', 'DisplayName', 'f_a (bvd)');
    fillmarkers; set(gca, 'ColorOrderIndex', 1);
    
    plot(x, [data.fr_peak]/1e9, 'o-', 'DisplayName', 'f_r (peak)'); 
    plot(x, [data.fa_peak]/1e9, 'o-', 'DisplayName', 'f_a (peak)');
    
    legend show;
    legend location NorthWest;
    xlabel 'reference k_t^2';
    ylabel 'frequency [GHz]';
end

function ktextplot(data)
    kt2ref = [data.kt2_ref];
    
    set(gca, 'ColorOrderIndex', 1);
    hb1 = ktplot(kt2ref, [data.kt2_bvd1],  'model 1 (peak)');
    hb2 = ktplot(kt2ref, [data.kt2_bvd2],  'model 2 (peak)');
    hb3 = ktplot(kt2ref, [data.kt2_bvd3],  'model 3 (peak)');
    fillmarkers;
    
    set(gca, 'ColorOrderIndex', 1);
    hp1 = ktplot(kt2ref, [data.kt2_peak1], 'model 1 (BvD)');
    hp2 = ktplot(kt2ref, [data.kt2_peak2], 'model 2 (BvD)');
    hp3 = ktplot(kt2ref, [data.kt2_peak3], 'model 3 (BvD)');
    
    plot(xlim, xlim, 'k--', 'DisplayName', 'Reference (x=y)'); 
    skiplegend;
    
    legend show; legend location SouthEast;

    xscale lin; xlabel 'reference k_t^2 [%]';
    yscale lin; ylabel 'extracted k_t^2 [%]';
end

function h = ktplot(kt2ref, kt2ext, name)
    h = plot(kt2ref*100, kt2ext*100, 'o-');
    h.DisplayName = name;
end

function data = hbar_compute(varargin)
    param = varargin;
    
    % main parameter
    name   = 'ePiezo';
    values = 1.4:0.2:2.4;
    
    f = [1e9:.1e6:5e9]';
    n = length(values);
    N = 400; % number of points 
    
    config = HBAR_loadconfig('HBAR_config_simple.txt');
    config = HBAR_parameters(config, param);
    
    for k = 1:n
        config = HBAR_parameters(config, {name, values(k)});
        [Zhbar, hbar]  = HBAR_v3(f, config);         
        [Zc, Zup, Zlo, F] = deripple(f, abs(Zhbar), N); 
        data(k) = ktextract(F, hbar, Zc, Zlo, Zup);
    end
end

function data = fbar_compute()
    % main parameter
    name   = 'ePiezo';
    values = 1.4:0.2:2.4;
    
    f = [5e9:.1e6:10e9]';
    n = length(values);
    N = 400; % number of points 
    
    fbar_param = {'tTopEl' 0 'tBotEl' 0 'tOxide' 0 'tSubst' 0};
    
    config = HBAR_loadconfig('HBAR_config_simple.txt');
    config = HBAR_parameters(config, fbar_param);
    
    for k = 1:n
        config = HBAR_parameters(config, {name, values(k)});
        [Zfbar, fbar]  = HBAR_v3(f, config);
        z = abs(Zfbar);
        data(k) = ktextract(f, fbar, z, z, z);
    end
end

function data = ktextract(F, hbar, Zc, Zlo, Zup)
    kt2_model1 = @(fr, fa) (fa^2 - fr^2)/fa^2;
    kt2_model2 = @(fr, fa) (pi^2/4)*(fr/fa)*(1-fr/fa);
    kt2_model3 = @(fr, fa) (pi/2*fr/fa)/tan(pi/2*fr/fa);
    
    % obtain fr, fa from peak detection
    %[fr_peak, fa_peak] = resonance_find(F, Zc);
    [fr_peak, ~] = resonance_find(F, Zlo);
    [~, fa_peak] = resonance_find(F, Zup);

    % obtain fr, fa from BvD fit 
    bvd    = bvdfit(F, Zc);
    fr_bvd = bvd.fr;
    fa_bvd = bvd.fa;
    
    data.hbar = hbar;
    data.bvd  = bvd;

    data.Zc = Zc;
    data.Zlo = Zlo;
    data.Zup = Zup;
    
    data.fr_peak = fr_peak;
    data.fa_peak = fa_peak;
    data.fr_bvd  = fr_bvd;
    data.fa_bvd  = fa_bvd;
    
    data.Zr_peak = interp1(F, Zlo, fr_peak);
    data.Za_peak = interp1(F, Zup, fa_peak);
    data.Zr_bvd  = interp1(F, bvd.zfit, fr_bvd);
    data.Za_bvd  = interp1(F, bvd.zfit, fa_bvd);

    data.kt2_ref = real(hbar.kt2);

    data.kt2_bvd1  = kt2_model1(fr_bvd, fa_bvd);
    data.kt2_bvd2  = kt2_model2(fr_bvd, fa_bvd);
    data.kt2_bvd3  = kt2_model3(fr_bvd, fa_bvd);

    data.kt2_peak1 = kt2_model1(fr_peak, fa_peak);
    data.kt2_peak2 = kt2_model2(fr_peak, fa_peak);
    data.kt2_peak3 = kt2_model3(fr_peak, fa_peak);
end
