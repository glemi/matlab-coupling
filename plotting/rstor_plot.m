function rstor_plot(f, Z)

    fig rstor:absz; clf;
    plot(f*1e-9, abs(Z))
    grid on
    title('Electrical impedance')
    xlabel('Frequency (GHz)')
    ylabel('|Z|')

    fig rstor:phasez; clf;
    plot(f*1e-9, rad2deg(angle(Z)), 'r')
    title('Phase')
    xlabel('Frequency (GHz)')
    ylabel('Angle')
    grid on

    fig rstor:smithGamma; clf;
    gamma = z2gamma(Z);
    smithchart(gamma)
    title('Smith Chart')
end