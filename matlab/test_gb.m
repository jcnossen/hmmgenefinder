% list proteins of e-coli
if exist('ecoli', 'var') == 0
    global ecoli;

    ecoli = GetGenbankFile('AE005174');
    disp 'Parsing features...';
    ecoli.fp = featuresparse(ecoli.Features);
end

numGenes = length(ecoli.fp.gene);
% set a subset of the genes here
%splitGene = 20;
splitGene = round(numGenes/2);
%endGene = 50;
endGene = numGenes;

clear testdata traindata;

% generate test set. The genome is split just after the gene with index 'splitGene'
disp 'generating test dataset';
testdata = GenomeCalc(ecoli, 'testdata');
range = [1 ecoli.fp.gene(splitGene).Indices(2)] ;
testdata.FilterDNA(range);
testdata.DispInfo();

% generate training set
disp 'generating training dataset';
traindata = GenomeCalc(ecoli, 'traindata');
range = [ecoli.fp.gene(splitGene).Indices(2)+1 ecoli.fp.gene(endGene).Indices(2)];
traindata.FilterDNA(range);
traindata.DispInfo();

genicLen = 0;
for i=1:length(ecoli.fp.gene)
    g = ecoli.fp.gene(i);
    genicLen = genicLen+max(g.Indices)-min(g.Indices);
end

disp(sprintf('Genic DNA length:%d. Total len: %d. Gene fraction: %f', ...
    genicLen, length(ecoli.Sequence), genicLen/length(ecoli.Sequence)));

testdata
traindata
