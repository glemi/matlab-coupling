function ripple_fit_analysis(data)
    fig ripplefit:eins; clf;
    global debug_plot; debug_plot = false;
    
    nostack = {'tSi', 1e-9, 'tSiO', 1e-9, 'tBotE', 1e-9, 'tTopE', 1e-9};
    [z, f] = rstor_model('Qpiezo', 10, 'fstep', 100e3, 'fstop', 4.5e9);
    %Zclean = rstor_model('fstep', 10e3, nostack{:});
   
    subplot(3,2,1);
    df = mean(diff(f(1:10)));
    ripple_distance = round(5e6/df);
    [Zu, Zl] = envelope(abs(z), ripple_distance, 'peak');
    Zc = sqrt(Zu.*Zl);
    
    subplot(3,2,2); cla;
    [fr, fa, fitdata] = ripplefit(f, abs(z));
    ok = [fitdata.valid];

    dfr = diff(fr);
    mdfr = median(dfr, 'omitnan');
    dfr(dfr > 1.01*mdfr | dfr < 0.98*mdfr) = NaN;
    
    fr1 = [fitdata.fr];
    dfr1 = diff(fr1);
    mdfr1 = median(dfr1, 'omitnan');
    dfr1(dfr1 > 1.01*mdfr1 | dfr1 < 0.98*mdfr1) = NaN;
   
    
    subplot(3,2,1);
    plot(fr(2:end)/1e9, dfr);
    plot(fr1(2:end)/1e9, dfr1);
    
    subplot(3,2,3);
    dualax left;
    plot(fr(ok)/1e9, [fitdata(ok).R], 'DisplayName', 'Damping R_m');
    plot(fr(ok)/1e9, [fitdata(ok).Q], 'DisplayName', 'Q');
    ylim([0 5000]);
    title 'R_m & Q';
    
    subplot(3,2,4);
    dualax left;
    plot(fr(ok)/1e9, [fitdata(ok).Cm], 'DisplayName', 'C_m');

    dualax right;
    plot(fr(ok)/1e9, [fitdata(ok).Lm], 'DisplayName', 'L_m');
    
    ylim([0 1e-3]);
    title 'C_m & L_m';
    %legend show; legend location best;
    
    subplot(3,2,5);
    plot(fr(ok)/1e9, [fitdata(ok).kt2]*100);
    
    title 'k_t^2 from BvD Fit';
end


% function [fr, fa, fit] = ripplefit(f, ir, ia, z)
%   
%     d = round(median(diff(ir)));
%     i = ir - round(d/3);
%     
%     n = length(i);
%     fr = zeros(1,n);  
%     dfr = zeros(1,n);
%     cr = zeros(1,n);
%    
%     irange = i(2):i(3);
%     fit(1) = bvdfit(f(irange), z(irange));
%     coeffs = fit(1).coeff_final;
% 
%     warning off stats:nlinfit:IterationLimitExceeded;
%     warning off stats:nlinfit:IllConditionedJacobian;
%     
%     for k = 2:1:n-1
%         irange = i(k):i(k+1);
%         frange = f(irange);
%         zrange = z(irange);
%             
%         fr(k) = NaN;
%         fa(k) = NaN;
%         fit(k) = bvdfit(frange, zrange, 'MaxIter', 30);
%         
%         if fit(k).valid
%             ffit = f(i(k)):20e3:f(i(k+1));
%             zfit = fit(k).model(ffit);
% 
%             [~, jr] = min(zfit);
%             [~, ja] = max(zfit);
%             fr(k) = ffit(jr);
%             fa(k) = ffit(ja);
%         
%             %% plot 
%             %continue;
%             cla;
%             plot(frange, zrange, '.');
%             ffit = f(i(k)):2e3:f(i(k+1));
%             zfit = fit(k).model(ffit);
%             plot(ffit, zfit);
%             plot(fr(k), fit(k).model(fr(k)), 'x');
%             plot(fa(k), fit(k).model(fa(k)), 'x');
%             ylim([min(zfit)*0.99 max(zfit)*1.01]);
%             %[min(zfit)*0.9 max(zfit)*1.1]
%             xlim([min(frange) max(frange)]);
%             set(gca, 'XTick', [], 'YTick', []);
%             title(sprintf('Fitting Ripples: %.0f%% (%.2fGHz)', 100*frange(1)./max(f), fr(k)/1e9));
%             drawnow;
%         end
%     end    
% %     fprintf('\n');
% %     fprintf('Q = %.2e\n', fit(k).Q);
% %     fprintf('Cm = %.2e\n', fit(k).Cm);
% %     fprintf('Lm = %.2e\n', fit(k).Lm);
% %     fprintf('C0 = %.2e\n', fit(k).C0);
% %     fprintf('C0 = %.2e\n', fit(k).Rm);
% 
%     warning on stats:nlinfit:IterationLimitExceeded;
%     warning on stats:nlinfit:IllConditionedJacobian;
% end
