% Calculates statistics for start and stop codons used in genes and their
% frequencies.
% Input:
%       seq - sequence with annotated genes
% Output:
%       res - statistics in form of a structure
function [res] = start_codon_stats_report(seq)
    n = length(seq.gene);
    for i = 1:length(HMM.Start_Codons)
        codon = HMM.Start_Codons{i};
        res.Start.(codon) = 0;
    end
    for i = 1:length(HMM.Stop_Codons)
        codon = HMM.Stop_Codons{i};
        res.Stop.(codon) = 0;
    end
    
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
                    res.Start.(startCodon) = res.Start.(startCodon) + 1;
                    res.Stop.(stopCodon) = res.Stop.(stopCodon) + 1;
                end
            end
        end
    end
    
    res.Start.Total = 0;
    for i = 1:length(HMM.Start_Codons)
        codon = HMM.Start_Codons{i};
        res.Start.Total = res.Start.Total + res.Start.(codon);
    end
    res.Stop.Total = 0;
    for i = 1:length(HMM.Stop_Codons)
        codon = HMM.Stop_Codons{i};
        res.Stop.Total = res.Stop.Total + res.Stop.(codon);
    end
    
    for i = 1:length(HMM.Start_Codons)
        codon = HMM.Start_Codons{i};
        res.Start.Percents.(codon) = res.Start.(codon) / res.Start.Total;
    end
    for i = 1:length(HMM.Stop_Codons)
        codon = HMM.Stop_Codons{i};
        res.Stop.Percents.(codon) = res.Stop.(codon) / res.Stop.Total;
    end
end