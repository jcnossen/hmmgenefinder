% Calculates statistics for codon usage in genes and their frequencies.
% Input:
%       seq - sequence with annotated genes
% Output:
%       res - statistics in form of a structure
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
function [str] = codon_stats_report_genic(seq)

    function [res] = sum_stats(a, b)
        fields = fieldnames(a);
        n = length(fields);
        for f = 1:n
            name = fields{f};
            res.(name) = a.(name) + b.(name);
        end
    end

    function [res] = output2csv()
        total = 0;
        
        endl = sprintf('\n');
        res = [';Train', endl];
        fields = fieldnames(codons);
        nFields = length(fields);
        
        for codonI = 1:nFields
            codonName = fields{codonI};
            total = total + codons.(codonName);
        end
        
        for codonI = 1:nFields
            codonName = fields{codonI};
            res = [res, codonName, ';', int2str(codons.(codonName)), ',', sprintf('%.4f', codons.(codonName) / total), ';', endl];
        end
        res = [res, 'Total;', int2str(total), ';'];
    end
    
    init = false;
    n = length(seq.gene);
    for i = 1:n
        st = seq.gene(i).Indices(1);
        fn = seq.gene(i).Indices(2);
        len = abs(fn - st) + 1;
        if (mod(len, 3) == 0)
            if (st < fn)
                gene = seq.Sequence(st:fn);
            else
                gene = seqrcomplement(seq.Sequence(fn:st));
            end
            
            startCodon = upper(gene(1:3));
            stopCodon = upper(gene(len - 2:len));
            if (HMM.is_start(startCodon) && HMM.is_stop(stopCodon))
                if (sum(nt2int(gene) > 4) == 0)
                    if (init)
                        codons = sum_stats(codons, codoncount(gene));
                    else
                        init = true;
                        codons = codoncount(gene);
                    end
                    codons.(startCodon) = codons.(startCodon) - 1;
                    codons.(stopCodon) = codons.(stopCodon) - 1;
                end
            end
        end
    end
    
    str = output2csv;
end