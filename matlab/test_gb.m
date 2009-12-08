% list proteins of e-coli
if exist('ecoli', 'var') == 0
    global ecoli;
    global ecoli_f; % parsed features

    ecoli = GetGenbankFile('AE005174');
    disp 'Parsing features...';
    ecoli_f = featuresparse(ecoli.Features);
end

%ListProteins(ecoli_f);

% ecoli_f
genicDNA = [];
intergenicDNA= [];
startIndex=1;
endIndex=length(ecoli_f.gene);
last=0;
lenSeq = length(ecoli.Sequence);
pos = 1;
for i=startIndex:endIndex
    gene = ecoli_f.gene(i);
    startPos=min(gene.Indices);
    endPos=max(gene.Indices);
    
    % add (pos:startPos) to intergenic
    if (startPos>pos)
        disp(sprintf('adding intergenic dna: start: %d, end: %d', pos,startPos));
        intergenicDNA = [intergenicDNA ecoli.Sequence(pos:startPos)];
        pos = endPos;
    end
    
    disp(sprintf('adding genic dna: start: %d, EndPos: %d', startPos, endPos));
    geneDNA = ecoli.Sequence(startPos:endPos);
    genicDNA = [genicDNA geneDNA];
    last=endPos;
end

if (pos<length(ecoli.Sequence))
    intergenicDNA = [intergenicDNA ecoli.Sequence(pos:lenSeq)];
end

disp(sprintf('Genic DNA length:%d. End: %d. Gene fraction: %f', length(genicDNA), last, length(genicDNA)/last));
disp(sprintf('Intergenic DNA length:%d.', length(intergenicDNA)));
