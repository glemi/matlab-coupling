% function pdfprintall(groupname)
% 
% Saves all figures in a group window as a pdf file. Uses a dialog
% to determine which directory to save the figurs. 
function pdfprintgroup(target, outdir)
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
    
    [fname, outdir] = uiputfile('*.pdf', 'Save Plot as PDF', outdir);
    pdffilename = fullfile(outdir, fname);
    
    if isnumeric(fname) || isempty(pdffilename)
        return;
    end
    
    n = length(figures);
    for k = 1:n
        hfigure = figures(k);
        fprintf('%s\n', hfigure.Name);
        filenames{k} = getfilename(hfigure, groupname);
        
        figure(hfigure);
        quickprint(filename, outdir);
    end
    
    try  
        delete(pdffilename);
        append_pdfs(pdffilename, filenames{:});
        
        for k = 1:n
            filepath = fullfile(filenames{k}, outdir);
            delete(filepath);
        end
    catch err
        errisp(err);
        message = sprintf('C');
        errordlg(err.message);
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


