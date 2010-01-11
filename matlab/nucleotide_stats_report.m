% Calculates statistics for nucleotide usage in given contig.
% Input:
%       seq - sequence with annotated genes
% Output:
%       res - statistics in form of a structure
function [res] = nucleotide_stats_report(seq)

    % Gets statistics for a single sequence
    function [res] = get_sequence_stats(seq)
        seqLen = length(seq);
        A = 0;
        C = 0;
        G = 0;
        T = 0;
        nuc = nt2int(seq);
        for j = 1:seqLen
            switch (nuc(j))
                case 1
                    A = A + 1;
                case 2
                    C = C + 1;
                case 3
                    G = G + 1;
                case 4
                    T = T + 1;
            end
        end
        Sum = A + C + G + T;
        Percents = struct('A', A / Sum, 'C', C / Sum, 'G', G / Sum, 'T', T / Sum);
        res = struct('A', A, 'C', C, 'G', G, 'T', T, 'Total', Sum, 'Percents', Percents);
    end

    % Sums up stats
    function [res] = sum_sequence_stats(a, b)
        A = a.A + b.A;
        C = a.C + b.C;
        G = a.G + b.G;
        T = a.T + b.T;
        Sum = a.Total + b.Total;
        Percents = struct('A', A / Sum, 'C', C / Sum, 'G', G / Sum, 'T', T / Sum);
        res = struct('A', A, 'C', C, 'G', G, 'T', T, 'Total', Sum, 'Percents', Percents);
    end

    res.Whole = get_sequence_stats(seq.Sequence);
    fprintf('[+] Got stats for whole genome.\n');
    res.Genic = get_sequence_stats('');
    res.Intergenic = get_sequence_stats('');

    geneC = 0;
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
                    res.Genic = sum_sequence_stats(res.Genic, get_sequence_stats(gene(4:len - 3)));
                    geneC = geneC + 1;
                end
            end
        end
    end
    fprintf('[+] Got stats for genic part: counted %i out of %i.\n', geneC, n);
    
    seq = HMM_Intergenic.get_intergenic(seq, false, false);
    fprintf('[i] Got intergenic sequences.\n');
    
    intC = 0;
    n = length(seq.gene);
    for i = 1:n
        st = seq.gene(i).Indices(1);
        fn = seq.gene(i).Indices(2);
        if (st < fn)
            gene = seq.Sequence(st:fn);
        else
            gene = seqrcomplement(seq.Sequence(fn:st));
        end
        if (sum(nt2int(gene) > 4) == 0)
            res.Intergenic = sum_sequence_stats(res.Intergenic, get_sequence_stats(gene));
            intC = intC + 1;
        end
    end
    fprintf('[+] Got stats for intergenic part: counted %i out of %i.\n', intC, n);
    fprintf('\n');
end