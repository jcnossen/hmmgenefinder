% list proteins of e-coli
if exist('ecoli', 'var') == 0
    global ecoli;
    global ecoli_f; % parsed features

    ecoli = GetGenbankFile('AE005174');
    disp 'Parsing features...';
    ecoli_f = featuresparse(ecoli.Features);
end

%ListProteins(ecoli_f);

% all features are stored in ecoli_f.gene
genicDNA = [];
startIndex=1;
endIndex=100;
last=0;
for i=startIndex:endIndex
    gene = ecoli_f.gene(i);
    startPos=min(gene.Indices);
    endPos=max(gene.Indices);

    disp(sprintf('Start: %d, EndPos: %d', startPos, endPos));
    geneDNA = ecoli.Sequence(gene.Indices(1):gene.Indices(2));
    genicDNA = [genicDNA geneDNA];
    last=endPos;
end
disp(sprintf('Genic DNA length:%d. End: %d. Gene fraction: %f', length(genicDNA), last, length(genicDNA)/last));
