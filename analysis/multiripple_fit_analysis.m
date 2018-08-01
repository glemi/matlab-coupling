function multiripple_fit_analysis(f, z, param, fitdata)
    
    global debug_plot; debug_plot = false;
    
    if nargin == 0
        f = (1e9:50e3:5e9)';
        %nostack = {'tSi', 1e-9, 'tSiO', 1e-9, 'tBotE', 1e-9, 'tTopE', 1e-9};
        [z, ~, param] = HBAR(f, 'tPiezo', 500e-9, 'tBotE', 200e-9, 'tTopE', 100e-9);
        
        % Ripple fit progress plot
        fitdata = multiRipplefit(f, abs(z), 3, 5e6, 'on');
    end
    
    
    fig ripplefit:data; clf;
    rawdataplot(f, z);
    
    
    fig ripplefit:results; clf;
    
    ok = [fitdata.valid];
    
    kt2ref = param.k2;
    cheeke = cheeke_compute(param);
    %Zclean = rstor_model('fstep', 10e3, nostack{:});
   
    df = mean(diff(f(1:10)));
    ripple_distance = round(5e6/df);
    [Zu, Zl] = envelope(abs(z), ripple_distance, 'peak');
    Zc = sqrt(Zu.*Zl);
    
    args = {'MinPeakDistance', ripple_distance};
    
    [~, ia] = findpeaks(abs(z), args{:});
    [~, ir] = findpeaks(1./abs(z), args{:});

    %findpeaks(abs(z11), args{:});
    
    fr = f(ir);
    fa = f(ia);
        
    %% delta-fr plot
    subplot(3,2,1);     
        dfr = diff(fr);
        mdfr = median(dfr, 'omitnan');
        dfr(dfr > 1.01*mdfr | dfr < 0.98*mdfr) = NaN;

        fr1 = [fitdata.fr];
        dfr1 = diff(fr1);
        mdfr1 = median(dfr1, 'omitnan');
        dfr1(dfr1 > 1.01*mdfr1 | dfr1 < 0.98*mdfr1) = NaN;

        dfr2 = [fitdata.dfr];
    
        plot(fr(1:end-1)/1e9, dfr, 'DisplayName', 'from Peak Detection');
        plot(fr1(1:end-1)/1e9, dfr1, 'DisplayName', 'diff of fr from fit');
        plot(fr1/1e9, dfr2, 'DisplayName', 'direct dfr from fit');
        
        %plot(xlim, [1 1]*cheeke.delta_fN, 'k-', 'LineWidth', 1); skiplegend;
        %plot(xlim, [1 1]*cheeke.delta_fT, 'k-', 'LineWidth', 1); skiplegend;
        
        legend show; legend location best;
        title '\Delta f_r resonance frequency spacing'
    
    %% Rm & Qm plot
    subplot(3,2,3);
        dualax left;
        plot(fr(ok)/1e9, [fitdata(ok).R], 'DisplayName', 'Damping R_m');
        plot(fr(ok)/1e9, [fitdata(ok).Q], 'DisplayName', 'Q');
        ylim([0 5000]);
        title 'R_m & Q';
    
    %% Cm & Lm plot
    subplot(3,2,4); title 'C_m & L_m';
        dualax left;
        plot(fr(ok)/1e9, [fitdata(ok).Cm], 'DisplayName', 'C_m');

        dualax right;
        plot(fr(ok)/1e9, [fitdata(ok).Lm], 'DisplayName', 'L_m');

        ylim([0 1e-3]);
        
        %legend show; legend location best;
    
    %% keff
    subplot(3,2,5); title 'k_{eff}^2 from BvD Fit [%]';
        keff = [fitdata(ok).kt2];
        plot(fr(ok)/1e9, keff*100);    
        
    
    %% keff & delta-fr 
    subplot(3,2,6);
        fr   = [fitdata.fr];
        keff = [fitdata.kt2];
        dfr  = [fitdata.dfr];
        
        [f0, dfr0, fit0result] = dfrvalleyfit(fr, dfr);
        [fm, dfrm, fitmresult] = dfrmaxfit(fr, dfr);
        keff0 = interp1(fr, keff, f0);

%         dualax left;
        dfr_analyze(fr, dfr, 'plot');
%         plot(fr, dfr); skipcolor;
%         plot(f0, dfr0, 'o'); fillmarkers;
%         axis manual;
%         plot(fit0result); 
%         plot(fitmresult); 
%         legend off;

%         dualax right;
%         plot(fr1, keff*100); skipcolor;
%         plot(f0, keff0*100, 'o'); fillmarkers; 
        
        kt2 = keff2kt(keff0);
        kt2 = cheeke.kt2Nfactor*keff0;
        [tS, dfrmax] = calcSubstThickness(fr, dfr);
        dfrvalleyfit(fr, dfr);
        
        fprintf('f0 = %.2f GHz\n', f0/1e9);
        fprintf('keff0 = %.2e\n', keff0);
        fprintf('kt2 = %.2f%% (actual kt2 = %.2f%%)\n', kt2*100, kt2ref*100);
        fprintf('tSi = %.1fum (dfr = %.2fMHz)\n', tS*1e6, dfrmax/1e6);
        
        
    subplot(3,2,2); cla;
        plot(fr, dfr); 
        plot(fr(2:end), dfr1); 
        
        
        plot(xlim, [1 1]*cheeke.delta_f0, 'k-', 'LineWidth', .5);
        plot(xlim, [1 1]*cheeke.delta_fN, 'k-', 'LineWidth', .5);
        plot(xlim, [1 1]*cheeke.delta_fT, 'k-', 'LineWidth', .5);
        
        plot([1 1]*cheeke.fN, ylim, 'k-', 'LineWidth', .5);
        plot([1 1]*cheeke.fT, ylim, 'k-', 'LineWidth', .5);
        
        
        x = fr1(1) + diff(xlim)*0.1; 
        dy = diff(ylim)*0.04;
        text(x, cheeke.delta_f0 + dy, '$\Delta f_0$', 'Interpreter', 'latex');
        text(x, cheeke.delta_fN + dy, '$\Delta f_N$', 'Interpreter', 'latex');
        
        
end


function [f0, dfr0, fitresult] = dfrvalleyfit(fr, dfr)
    width = 0.1e9;
    [~, i] = min(dfr);
    f1 = fr(i) - width/2;
    f2 = fr(i) + width/2;
    i = (fr > f1 & fr < f2);
    %  p1*x^2 + p2*x + p3
    fitresult = fit(fr(i)', dfr(i)', 'poly2');
    a = fitresult.p1;
    b = fitresult.p2;
    
    f0 = -b/(2*a);
    dfr0 = feval(fitresult, f0);
end

function [f0, dfr0, fitresult] = dfrmaxfit(fr, dfr)
    width = 0.2e9;
    [~, i] = min(dfr);
    [~, i] = max(dfr(1:i));
    f1 = fr(i) - width/2;
    f2 = fr(i) + width/2;
    i = (fr > f1 & fr < f2);
    %  p1*x^2 + p2*x + p3
    fitresult = fit(fr(i)', dfr(i)', 'poly2');
    a = fitresult.p1;
    b = fitresult.p2;
    
    f0 = -b/(2*a);
    dfr0 = feval(fitresult, f0);
end

function [tS, dfr0] = calcSubstThickness(fr, dfr)    
    [~, dfr0] = dfrmaxfit(fr, dfr);

    vS = 8.44065e+03;
    tS = vS/dfr0/2;
end

function kt2 = keff2kt(keff2)
    tPiezo  = 500e-9;  % AlScN layer thickness
    tTopE   = 100e-9;   % Top electrode thickness
    tBotE   = 200e-9;  % Bottom electrode thickness
    tSiO2   = 200e-9;  % SiO2 thickness
    tSi     = 725e-6;  % Si thickness
    
    denPiezo = 3280; % AlScN mass density
    denSiO2  = 2197;  % Mass density of SiO2
    denSi   = 2330;  % Mass density of Si
    denPt   = 21450; % Mass density of electrodes (Pt)
    
    mPiezo = tPiezo*denPiezo;
    mTopE  = tTopE*denPt;
    mBotE  = tBotE*denPt;
    mSi    = tSi*denSi;
    mSiO2  = tSiO2*denSiO2;
    mSubst = mSi+mSiO2;
    
    vSubst = 8.4407e+03;
    vPiezo = 9.6583e+03;
    
    R = vPiezo/vSubst * (tSi + .5*mBotE/denSi)/(tPiezo + mTopE/denPiezo + .5*mTopE/denPiezo);
    
    x = mPiezo *( mPiezo + mSubst + mTopE + mBotE ) / (mPiezo + mTopE + 1/2*mBotE)^2;
    kt2 = keff2*x;
end
