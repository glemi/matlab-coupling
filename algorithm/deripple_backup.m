function [z, zupper, zlower] = deripple(f, z, fi)
    %global debug_plot;

    if ~isreal(z)
        z = abs(z);
    end
    
    [~, i] = findpeaks(z);
    [~, j] = findpeaks(1./z);

    if nargin < 3
        fi = f;
    end
    
    zupper = interp1(f(i), z(i), fi);
    zlower = interp1(f(j), z(j), fi);
    zcentr = sqrt(zupper.*zlower);

    z = zcentr;
    
%     if debug_plot
%         plot(f1, zcentr, 'w', 'LineWidth', 1);
%         plot(f1, zupper);
%         plot(f1, zlower);
% 
%         [zr, ir] = min(zcentr);
%         [za, ia] = max(zcentr);
%         
%         yscale log;
%     end
end