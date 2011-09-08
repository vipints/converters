function reads = mat_cell2struct(filename)
    try
        l = load(filename);
    catch 
		fprintf('could not load pre-processed reads information in .mat format %s\n', filename)
        quit; 
    end
    x = 1;
    for i=1:length(l.read_dist) 
        read = l.read_dist{x};  
        if ~isempty(read.Oth),
            read.Oth = double(read.Oth);
            read.Oth = reshape(read.Oth, 1, length(read.Oth));
            read.Oth_MM = double(read.Oth_MM);
            read.Oth_MM = reshape(read.Oth_MM, 1, length(read.Oth_MM));
        end;
        if ~isempty(read.TE),
            read.TE = double(read.TE);
            read.TE = reshape(read.TE, 1, length(read.TE));
            read.TE_MM = double(read.TE_MM);
            read.TE_MM = reshape(read.TE_MM, 1, length(read.TE_MM));
        end;
        if ~isempty(read.ribo),
            read.ribo = double(read.ribo);
            read.ribo = reshape(read.ribo, 1, length(read.ribo));
            read.ribo_MM = double(read.ribo_MM);
            read.ribo_MM = reshape(read.ribo_MM, 1, length(read.ribo_MM));
        end;
        if ~isempty(read.CG),
            read.CG = double(read.CG);
            read.CG = reshape(read.CG, 1, length(read.CG));
            read.CG_MM = double(read.CG_MM);
            read.CG_MM = reshape(read.CG_MM, 1, length(read.CG_MM));
        end;
        read.count = double(read.count);
        read.orient = reshape(read.orient, 1, length(read.orient));
		reads(length(l.read_dist)-i+1) = read; 
        x = x + 1;
    end
    save('-v7', filename, 'reads');
    clear l
