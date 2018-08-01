function h = latexfigure(figname, text)
    hf = fig(figname); clf;
    hf.MenuBar = 'none';
    hf.ToolBar = 'none';
    axis off;
    
    label([.5 .5], text, 'HorizontalAlignment', 'center');

end

