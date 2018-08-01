function [data] = read_s1p(filename, title)
    % http://cp.literature.agilent.com/litweb/pdf/ads2004a/cktsim/ck04a12.html
    % http://cp.literature.agilent.com/litweb/pdf/ads2004a/cktsim/ck04a8.html

    % Format string for each line of text:
    %formatSpec = '%f%f%f%f%f%f%f%f%f%[^\n\r]';
    formatSpec = '%f%f%f';
    
    % Open the text file.
    fileID = fopen(filename,'r');

    %% Header
    cont = true;
    iLine = 0;
    while cont
        iLine = iLine + 1;
        line = fgetl(fileID);
        
        if strncmp(line, '! VAR ', 6)
            header = strsplit(line(7:end), '=');
            name = lower(header{1});
            value = header{2};
            num = str2double(value);
            if isnan(num)
                data.(name) = value;
            else
                data.(name) = num;
            end
        elseif strncmp(line, 'VAR ', 4)
            header = strsplit(line(5:end), '=');
            name = lower(header{1});
            value = header{2};
            num = str2double(value);
            if isnan(num)
                data.(name) = value;
            else
                data.(name) = num;
            end
        elseif line(1) == '!'
            % do nothing
        elseif strcmp(line(1:5), 'BEGIN')
            % do nothing
        elseif line(1) == '#'
            strings = strsplit(line(3:end));
            data.freq_unit = strings{1};
            data.signal_parmType = strings{2};
            data.signal_parmFormat = strings{3};
            data.sys_impedance = str2double(strings{5});
        elseif ~isempty(str2num(line))
            nHeaderLines = iLine-1;
            cont = false;
        end
    end
    
    %% filename
    [~, name, ext] = fileparts(filename);
    data.filename = [name ext];
    
    %% title
    if nargin >= 2
        data.title = title;
    end
    
    %% Read columns of data according to format string.
    frewind(fileID);
    dataArray = textscan(fileID, formatSpec, Inf, 'HeaderLines', nHeaderLines, 'ReturnOnError', false);

    %% Close the text file.
    fclose(fileID);

    %% Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post
    % processing code is included. To generate code which works for
    % unimportable data, select unimportable cells in a file and regenerate the
    % script.

    %% Allocate imported array to column variable names
    
    % Raw Data:
    F     = dataArray{:, 1};
    N11_1 = dataArray{:, 2};
    N11_2 = dataArray{:, 3};
    
    data.F     = F;
    data.N11_1 = N11_1;
    data.N11_2 = N11_2;
    
    switch data.freq_unit(1:2)
        case 'Hz', data.f = dataArray{:, 1};
        case 'KH', data.f = dataArray{:, 1} * 1e3;
        case 'MH', data.f = dataArray{:, 1} * 1e6;
        case 'GH', data.f = dataArray{:, 1} * 1e9;
        otherwise, data.f = dataArray{:, 1};
    end

    switch data.signal_parmFormat
        case 'RI' % real / imaginary
            data.s11 = complex( N11_1, N11_2);
        case 'MA' % magnitude / angle (degrees!)
            data.s11 = N11_1.*exp(1j*N11_2*pi/180);
        case 'DB' % decibels / angle (degrees!)
            data.s11 =  db2mag(N11_1).*exp(1j*N11_2*pi/180); %10.^(N11_1/20).*exp(1j*N11_2);
    end

    data.n = length(data.f);
    data.s = zeros(2,2,data.n);
    data.s(1, 1, :) = data.s11;
    data.z = s2z(data.s, data.sys_impedance);
    data.y = s2y(data.s, data.sys_impedance);
end