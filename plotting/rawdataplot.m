function rawdataplot(f, Z)
    fGHz = f/1e9;

    subplot(2,2,1);
    plot(fGHz, abs(Z));
    xscale log; yscale log;
    title '|Z| [\Omega]';
    xlabel 'frequency [GHz]';
    %ylim([5 150]);
    
    subplot(2,2,2);
    plot(fGHz, angle(Z)*180/pi);
    title '\phi(Z) [\circ]';
    xlabel 'frequency [GHz]';
    %ax.YLim(1) = -90;
    
    ax = subplot(2,2,3);
    plot(fGHz, real(Z));
    title 'Re\{Z\}';
    xlabel 'frequency [GHz]';
    %ax.YLim(1) = 0;
    
    subplot(2,2,4);
    plot(fGHz, -imag(Z));
    xscale log; yscale log;
    title '-Im\{Z\}';
    xlabel 'frequency [GHz]';
    %ylim([0 120]);
end