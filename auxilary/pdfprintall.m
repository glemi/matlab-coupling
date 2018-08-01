% function pdfprintall(groupname)
% 
% Saves all figures in a group window as a pdf file. Uses a dialog
% to determine which directory to save the figurs. 
function pdfprintall(target, outdir)
    opt groupname char '';
    opt outdir char 'figures';
    
    if isempty(target)
        target = gcf;
    end
    
    if ischar(target)
        groupname = target;
    elseif ishandle(target)
        figname = get(target, 'Name');
        parts = strsplit(figname, ':');
        if length(parts) == 2
            groupname = parts{1};
        else
            return;
        end
    else
        error 'invalid argument';
    end
        
    figures = findobj('Type', 'figure');
    fignames = {figures.Name};
    i = strncmp([groupname ':'], fignames, length(groupname)+1);
    figures = figures(i);
    
    outdir = getdir(outdir);
    if isnumeric(outdir)
        return;
    end
    
    n = length(figures);
    for k = 1:n
        hfigure = figures(k);
        fprintf('%s\n', hfigure.Name);
        filename = getfilename(hfigure, groupname);
        
        figure(hfigure);
        quickprint(filename, outdir);
    end
    
end

function filename = getfilename(hfigure, groupname)
    filename = hfigure.Name;
    if ~isempty(groupname)
        %filename = filename(length(groupname)+2:end);
        filename = strrep(filename, ':', '.');
    end
    filename = [filename '.pdf'];
end

function outdir = getdir(startdir)
    if isempty(startdir) || ~isdir(startdir)
        if isdir('figures')
            startdir = fullfile(pwd, 'figures');
        else
            startdir = pwd;
        end
    end
    outdir = uigetdir(startdir, 'pick directory');
end


