% Returns sequences for training HMM (transforms annotated genes
% into cell array of sequences)
% Input:
%       seq - sequence with annotated genes
% Output:
%       res - cell array of sequences (to be fed to Viterbi
%             algorithm).
function [res] = get_training_sequences(seq)
    n = length(seq.gene);
    res = cell(0);
    resCount = 0;
    for i = 1:n
        s = min(seq.gene(i).Indices);
        f = max(seq.gene(i).Indices);
        str = seq.Sequence(s:f);
        if (seq.gene(i).Indices(1) > seq.gene(i).Indices(2))
            str = seqrcomplement(str);
        end
        %int_str = nt2int(str);
        %if (sum(int_str > 4) == 0) % Traditional only
        resCount = resCount + 1;
        res{resCount} = num2cell(str);
        %end
    end
end
