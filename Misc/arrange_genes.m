function genes = arrange_genes(filename)
    try
        l = load(filename);
    catch 
		fprintf('could not load pre-processed annotation in .mat format %s\n', filename)
        quit; 
    end
    x = 1;
    for i=1:length(l.gene_models) 
        gene = l.gene_models{x};  
        % arranging transcript_id
        gene.transcripts = reshape(gene.transcripts, 1, length(gene.transcripts));
        % reshaping exons
        for j=1:length(gene.exons)
            gene.exons{j} = double(gene.exons{j});
        end    
        gene.exons = reshape(gene.exons, 1, length(gene.exons));
        gene.id = double(gene.id);
		genes(length(l.gene_models)-i+1) = gene; 
        x = x + 1;
    end
    save('-v7.3', filename, 'genes');
    clear l
