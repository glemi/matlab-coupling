function [hline, htext] = vline(x, y, labeltext)
    
    yl = ylim; xl = xlim;

    hline = plot([x x], yl, 'k-');
    skiplegend;
    
    dx = diff(xl)*0.01;
    y  = diff(yl)*y + yl(1);
    htext = text(x+dx, y, labeltext);
    
    %htext.BackgroundColor = 'w';
    if length(strfind(labeltext, '$')) >= 2
        htext.Interpreter = 'latex';
    end
    
end