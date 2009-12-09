classdef GenomeCalc < handle
    properties
        range
        intergenic
        genic
        genome
        name
    end
    
    methods
        % constructor
        function gd = GenomeCalc(gbk, name)
            gd.genome = gbk;
            gd.name = name;
        end
        % Splits DNA in a genic part and an intergenic part
        % Parameters:
        % genome = the GBK object as returned by GetGenbankFile
        % featureData = the parsed feature data
        % range = range of dna to filter in array, 
        %         like: [startNucleotide endNucleotide]
        % 
        % NOTE: this function does not split genes in 2. If the range ends before a
        % gene ends, it will copy beyond the range and include the full gene
        function FilterDNA(inst, range)
            genic = [];
            intergenic = [];
            
            % parsed features struct (see featuresparse function)
            fp = inst.genome.fp;
            inst.range = range;

            % decide at which gene to start and stop
            if isempty(range)
                startPos = 1;
                endPos = length(inst.genome.Sequence);
            else
                startPos = range(1);
                endPos = range(2);
            end


            % find the first gene within the given range, and store its index in geneIndex
            geneIndex=1;
            for g=1:length(fp.gene)
                gene = fp.gene(g);

                % copy gene location
                % complementary strand has reversed indices, so we can't just use Indices(1)
                gloc = [min(gene.Indices) max(gene.Indices)]; %
                if (gloc(1) >= startPos)
                    geneIndex = g;
                    break
                end
            end

            % loop through all genes
            pos = startPos; % nucleotide pos
            while pos <= endPos
                gene = fp.gene(geneIndex);
                gloc = [min(gene.Indices) max(gene.Indices)];

                if gloc(2) > endPos
                    break
                end

                % add (pos:startPos) to intergenic
                if (gloc(1)>pos)
                    %disp(sprintf('adding intergenic dna: start: %d, end: %d', pos,gloc(1)-1));
                    intergenic = [intergenic inst.genome.Sequence(pos:gloc(1)-1)];
                end

                %disp(sprintf('adding genic dna: start: %d, EndPos: %d', gloc(1), gloc(2)));
                geneDNA = inst.genome.Sequence(gloc(1):gloc(2));
                genic = [genic geneDNA];
                pos = gloc(2)+1; % set new pos to just after the gene
                geneIndex=geneIndex+1;
            end

            % add the piece of intergenic DNA that follows the last gene
            if (pos<=endPos)
                %disp(sprintf('adding intergenic dna: start: %d, end: %d', pos,endPos));
                intergenic = [intergenic inst.genome.Sequence(pos:endPos)];
            end
            
            inst.genic = genic;
            inst.intergenic = intergenic;
        end
    
        function DispInfo(gc)
            len = gc.range(2)-gc.range(1);
            disp(sprintf('%s: total length: %f. Genic part: %f. Intergenic part: %f', ...
                gc.name, len, length(gc.genic)/len, length(gc.intergenic)/len));
        end
        
    end

    
        
end

