addpath ~/svn/projects/mGene_core/data_processing_signals/ 
addpath ~/svn/projects//mGene_core/auxiliary_data/
%addpath /fml/ag-raetsch/home/jonas/svn/projects/mGene_core/shared_tools/ngs/

fn_bams={'s2_replicate_mmr_uniq_sort.bam',  's3_replicate_mmr_uniq_sort.bam',  's5_replicate_mmr_uniq_sort.bam',  's6_replicate_mmr_uniq_sort.bam', ...
         'Lane2_mmr_uniq_sort.bam', 'Lane3_mmr_uniq_sort.bam',  'Lane4_mmr_uniq_sort.bam'};

descr_bam={'negative control (WT w/ untagged Scc1)', 'positive control (WT w/ tagged Scc1)', 'Scc3 mutation strain',  'PDS5 mutation strain R1', ...
           'PDS5 mutation strain R2', 'PDS5 mutation strain R3', 'deletion of rad61 gene', 'PDS5 mutation strain',  'masked positive control/PDS5', 'masked positive control/PDS5 over negative control', 'A', 'B'} ;

chrs={'chrI', 230208;
'chrII',   813178;
'chrIII',  316617;
'chrIV',   1531919;
'chrV',    576869;
'chrVI',   270148;
'chrVII',  1090947;
'chrVIII', 562643;
'chrIX',   439885;
'chrX',    745741;
'chrXI',   666454;
'chrXII',  1078175;
'chrXIII', 924429;
'chrXIV',  784334;
'chrXV',   1091289;
'chrXVI',  948062;
'chrMito', 85779;
'2-micron', 6318};
%chrs={'chrIII', 150000};

fd=fopen('/fml/ag-raetsch/nobackup2/projects/sequencing_runs/yeast/reads/saccharomyces_cerevisiae_gene.gff', 'r') ;
%fd=fopen('chrIII.bed', 'r') ;
genes{18,2}=[] ;
%genes{1,2}=[] ;

while ~feof(fd),
    line=fgetl(fd) ;
    if ~ischar(line), break ; end ;
    items = separate(line) ;
    chr_idx=strmatch(items{1}, chrs(:,1), 'exact') ;
    if isequal(items{4}, '+'),
        strand = 1 ;
    elseif  isequal(items{4}, '-'),
        strand = 2 ;
    else
        error('strand') ;
    end ;
    genes{chr_idx, strand}=[genes{chr_idx, strand}; str2num(items{2}) str2num(items{3})] ;
end ;
fclose(fd) ;

for k=1:size(chrs,1),

    k

    gg=[] ;
    gg.chr=chrs{k,1} ;
    gg.start=1;
    gg.stop=chrs{k,2} ;
    gg.strand='+'; 

    window=500 ;

    maxval=inf ;
    conf_filter=[] ;
    conf_filter.mismatch=1 ;
    conf_filter.exon_len=1 ;
    conf_filter.intron=1 ;

    tracks=[] ;
    for f=1:length(fn_bams),
        gg.tracks=[] ;
        gg = add_reads_from_bam(gg, fn_bams{f}, 'exon_track', '', maxval, conf_filter);
        track=gg.tracks ;
    
        w_idx=window/2+1:length(track)-window/2 ;
        w_idx=w_idx(1:100:end) ;
        track=zeros(1,length(w_idx)) ;
        for i=1:length(w_idx),
            if mod(i,1000)==0, i, end ;
            idx=(w_idx(i)-window/2):(w_idx(i)+window/2) ;
            track(i)=mean(gg.tracks(idx)) ;
        end ;

        tracks=[tracks; track/median(track)] ;
    end ;

    logtracks=[] ;
    for i=2:size(tracks,1),
        logtracks(i,:)=log(tracks(i,:)./tracks(1,:)) ;
    end ;

    figure(1) ;
    clf
    subplot(2,1,1) ;
    logtracks(8,:)=log(((tracks(5,:)+tracks(6,:))/2)./tracks(1,:)) ;
    plot_idx=[2 8] ;
    plot(w_idx, logtracks(plot_idx,:)','LineWidth',1) ;
    legend(descr_bam(plot_idx)) ;
    x=axis ;
    x(2)=chrs{k,2} ;
    x(3)=-2 ;
    x(4)=2 ;
    ylabel('log-fold change') ;
    xlabel('genome position') ;
    grid on
    axis(x) ;

    hold on ;
    for q=1:size(genes{k,1},1)
        plot(genes{k,1}(q, :), [-1.8 -1.8], 'b-','LineWidth',1) ;
        plot(genes{k,1}(q, 1)*[1 1], [-1.8 0], 'b:','LineWidth',0.1) ;
        plot(genes{k,1}(q, 2)*[1 1], [-1.8 0], 'b:','LineWidth',0.1) ;
    end;
    for q=1:size(genes{k,2},1)
        plot(genes{k,2}(q, :), [-1.85 -1.85], 'g-','LineWidth',1) ;
        plot(genes{k,2}(q, 1)*[1 1], [-1.85 0], 'g:','LineWidth',0.1) ;
        plot(genes{k,2}(q, 2)*[1 1], [-1.85 0], 'g:','LineWidth',0.1) ;
    end;

    title(sprintf('%s (contig %i)', chrs{k,1} ,k)) ;

    subplot(2,1,2) ;
    logtracks(9,:)=-log(((tracks(5,:)+tracks(6,:))/2)./tracks(2,:)) ;
    %logtracks(9,logtracks(2,:)>0.5 & logtracks(8,:)>0.5)=0 ;
    logtracks(9,logtracks(2,:)<0.25 & logtracks(8,:)<0.25)=0 ;
    logtracks(10,logtracks(2,:)>0.5 & logtracks(8,:)>0.5)=0 ;
    logtracks(10,:)=-log(((tracks(5,:)+tracks(6,:))/2)./tracks(2,:)) ;
    logtracks(10,logtracks(2,:)<0.25 & logtracks(8,:)<0.25)=0 ;
    logtracks(10,logtracks(2,:)<0.25 & logtracks(8,:)>0.25)=0 ;
    logtracks(10,logtracks(2,:)>0.25 & logtracks(8,:)>0.25)=0 ;
    logtracks(11,:)=log(((tracks(5,:)+tracks(6,:))/2)./tracks(2,:)) ;
    logtracks(11,logtracks(2,:)<0.25 & logtracks(8,:)<0.25)=0 ;
    logtracks(11,logtracks(2,:)>0.25 & logtracks(8,:)<0.25)=0 ;
    logtracks(11,logtracks(2,:)>0.25 & logtracks(8,:)>0.25)=0 ;
    %plot_idx=[2 8 9 10 11] ;
    plot_idx=[2 8 9] ;
    plot(w_idx, logtracks(plot_idx,:)','LineWidth',1) ;
    legend(descr_bam(plot_idx)) ;
    x=axis ;
    x(2)=chrs{k,2} ;
    x(3)=-2 ;
    x(4)=2 ;
    ylabel('log-fold change') ;
    xlabel('genome position') ;
    grid on
    axis(x) ;

    hold on ;
    for q=1:size(genes{k,1},1)
        plot(genes{k,1}(q, :), [-1.8 -1.8], 'b-','LineWidth',1) ;
        plot(genes{k,1}(q, 1)*[1 1], [-1.8 0], 'b:','LineWidth',0.1) ;
        plot(genes{k,1}(q, 2)*[1 1], [-1.8 0], 'b:','LineWidth',0.1) ;
    end;
    for q=1:size(genes{k,2},1)
        plot(genes{k,2}(q, :), [-1.85 -1.85], 'g-','LineWidth',1) ;
        plot(genes{k,2}(q, 1)*[1 1], [-1.85 0], 'g:','LineWidth',0.1) ;
        plot(genes{k,2}(q, 2)*[1 1], [-1.85 0], 'g:','LineWidth',0.1) ;
    end;

    if 0,
        figure(1) ;
        subplot(3,1,3) ;
        logtracks(10,:)=log(tracks(2,:)./((tracks(5,:)+tracks(6,:))/2)) ;
        logtracks(10,logtracks(2,:)<0.25 & logtracks(8,:)<0.25)=0 ;
        plot_idx=[2 8 10] ;
        plot(w_idx, logtracks(plot_idx,:)') ;
        legend(descr_bam(plot_idx)) ;
        x=axis ;
        x(2)=chrs{k,2} ;
        x(3)=-2 ;
        x(4)=2 ;
        grid on
        axis(x) ;
    end ;

    fn_out=sprintf('%s.eps', chrs{k,1}) 
    print(fn_out, '-depsc2')
end ;
