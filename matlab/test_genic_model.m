% Runs viterbi algorithm for genes and reports model performance (number of
% correctly modeled genes)
% Input:
%       <hmm> - HMM object to test
%       <seq> - sequence with annotated genes
% Output:
%       [res] - testing results in form of a structure
function [res] = test_genic_model(hmm, seq)
    hmm = HMM(hmm);
    hmm.insert_state(1, 'entry_point', [0 0 0 0]);
    hmm.add_edge(1, 'start_codons_AGT', 1);
    
    n = length(seq.gene);
    seq.Sequence = upper(seq.Sequence);
    
    GoodGenes = 0;
    CorrectlyFound = 0;
    CorrectlyFoundLen3 = 0;
    CorrectlyFoundNotLen3 = 0;
    Len3 = 0;
    NotLen3 = 0;
    for i = 1:n
        st  = seq.gene(i).Indices(1);
        fn  = seq.gene(i).Indices(2);
        len = abs(fn - st) + 1;
        if (st < fn)
            gene = seq.Sequence(st:fn);
        else
            gene = seqrcomplement(seq.Sequence(fn:st));
        end
        startCodon = gene(1:3);
        stopCodon = gene(len - 2:len);
        if (HMM.is_start(startCodon) && HMM.is_stop(stopCodon))
            if (sum(nt2int(gene) > 4) == 0)
                GoodGenes = GoodGenes + 1;
                states = hmm.viterbi(gene);
                good = (sum(strncmp('start', states(1:3), 5)) == 3) && (sum(strncmp('stop', states(len - 2:len), 4)) == 3);
                if (good)
                    CorrectlyFound = CorrectlyFound + 1;
                end
                if (mod(len, 3) == 0)
                    if (good)
                        CorrectlyFoundLen3 = CorrectlyFoundLen3 + 1;
                    end
                    Len3 = Len3 + 1;
                else
                    if (good)
                        CorrectlyFoundNotLen3 = CorrectlyFoundNotLen3 + 1;
                    end
                    NotLen3 = NotLen3 + 1;
                end
            end
        end
    end
    Percents = struct('CorrectlyFound', CorrectlyFound / GoodGenes, 'CorrectlyFoundLen3', CorrectlyFoundLen3 / Len3, 'CorrectlyFoundNotLen3', CorrectlyFoundNotLen3 / NotLen3);
    res = struct('TotalGenes', n, 'GoodGenes', GoodGenes, 'Len3', Len3, 'NotLen3', NotLen3, 'CorrectlyFound', CorrectlyFound, 'CorrectlyFoundLen3', CorrectlyFoundLen3, 'CorrectlyFoundNotLen3', CorrectlyFoundNotLen3, 'Percents', Percents);
end