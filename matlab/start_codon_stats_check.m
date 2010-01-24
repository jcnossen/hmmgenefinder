% Counts the number of start and stop codons that occur in the middle of a
% gene.
% Input:
%       <seq> - sequence with annotated genes
% Output:
%       <res> - statistics in form of a structure
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
function [res] = start_codon_stats_check(seq)

    function [res] = sum_stats(a, b)
        for iR = 1:length(HMM.Start_Codons)
            codon = HMM.Start_Codons{iR};
            res.Start.(codon) = a.Start.(codon) + b.(codon);
        end
        for iR = 1:length(HMM.Stop_Codons)
            codon = HMM.Stop_Codons{iR};
            res.Stop.(codon) = a.Stop.(codon) + b.(codon);
        end
    end
    
    for i = 1:length(HMM.Start_Codons)
        codon = HMM.Start_Codons{i};
        res.Start.(codon) = 0;
    end
    for i = 1:length(HMM.Stop_Codons)
        codon = HMM.Stop_Codons{i};
        res.Stop.(codon) = 0;
    end

    n = length(seq.gene);
    geneC = 0;
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
                    geneC = geneC + 1;
                    res = sum_stats(res, codoncount(gene));
                    res.Start.(startCodon) = res.Start.(startCodon) - 1;
                    res.Stop.(stopCodon) = res.Stop.(stopCodon) - 1;
                end
            end
        end
    end
    fprintf('Processed %i / %i genes.\n', geneC, n);
end