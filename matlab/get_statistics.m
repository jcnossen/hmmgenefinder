% Returns statistics of nucleotide and codon usage in given DNA regions.
% Used by HMM family classes when constructing models.
% Input:
%      <seq>                    - sequence with annotated regions
%      <traditionalNucleotides> - boolean switch to gather statistics only
%                                 for genes with traditional nucleotides
%      <length3>                - boolean switch to gather statistics only
%                                 for genes with length devisible by 3
%      <traditionalCodons>      - boolean switch to gather statistics only
%                                 for genes with traditional start and stop
%                                 codons
%      <countStart>             - boolean switch to include start and stop
%                                 codon statistics in results (if false,
%                                 nucleotide and codon statistics on start
%                                 and stop codons is not counted).
% Output:
%      <sum_arr>                - array with codon statistics (as returned
%                                 by codoncount)
%      <sum_nuc>                - array with nucleotide statistics
%      <geneC>                  - number of regions processed (counts only
%                                 regions answering requested contraints)
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
function [sum_arr, sum_nuc, geneC] = get_statistics(seq, traditionalNucleotides, length3, traditionalCodons, countStart)
    if (nargin < 2)
        traditionalNucleotides = false;
    end

    n = length(seq.gene);
    sum_arr = zeros([4 4 4]);
    sum_nuc = zeros(1, 16);
    geneC = 0;
    for i = 1:n
        s = min(seq.gene(i).Indices);
        f = max(seq.gene(i).Indices);
        len = f - s + 1;
        if (~length3 || mod(len, 3) == 0)
            str = seq.Sequence(s:f);
            if (seq.gene(i).Indices(1) > seq.gene(i).Indices(2))
                str = seqrcomplement(str);
            end
            if (len >= 3)
                startCodon = upper(str(1:3));
                stopCodon = upper(str(len - 2:len));
            else
                startCodon = '';
                stopCodon = '';
            end
            if (~traditionalCodons || (HMM.is_start(startCodon) && HMM.is_stop(stopCodon)))
                int_str = nt2int(str);
                if (~traditionalNucleotides || sum(int_str > 4) == 0)
                    geneC = geneC + 1;
                    [tmp, tmp_arr] = codoncount(str);
                    sum_arr = sum_arr + tmp_arr;
                    len = length(str);
                    for j = 1:len
                        ind = int_str(j);
                        sum_nuc(ind) = sum_nuc(ind) + 1;
                    end
                    if (~countStart)
                        if (length(startCodon) == 3)
                            sum_arr(nt2int(startCodon(1)), nt2int(startCodon(2)), nt2int(startCodon(3))) = sum_arr(nt2int(startCodon(1)), nt2int(startCodon(2)), nt2int(startCodon(3))) - 1;
                        end
                        if (length(stopCodon) == 3)
                            sum_arr(nt2int(stopCodon(1)), nt2int(stopCodon(2)), nt2int(stopCodon(3))) = sum_arr(nt2int(stopCodon(1)), nt2int(stopCodon(2)), nt2int(stopCodon(3))) - 1;
                        end
                        sum_nuc(nt2int(startCodon(1))) = sum_nuc(nt2int(startCodon(1))) - 1;
                        sum_nuc(nt2int(startCodon(2))) = sum_nuc(nt2int(startCodon(2))) - 1;
                        sum_nuc(nt2int(startCodon(3))) = sum_nuc(nt2int(startCodon(3))) - 1;
                        sum_nuc(nt2int(stopCodon(1))) = sum_nuc(nt2int(stopCodon(1))) - 1;
                        sum_nuc(nt2int(stopCodon(2))) = sum_nuc(nt2int(stopCodon(2))) - 1;
                        sum_nuc(nt2int(stopCodon(3))) = sum_nuc(nt2int(stopCodon(3))) - 1;
                    end
                end
            end
        end
    end
end