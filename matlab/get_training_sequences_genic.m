% Returns sequences for training HMM (used to export data for GHMM)
% Input:
%       seq - sequence with annotated genes
% Output:
%       res - cell array with genic regions
function [res] = get_training_sequences_genic(seq)
    % Just in case, but better do it when loading the genome
    seq.Sequence = upper(seq.Sequence);
    
    seq = HMM.select_genes(seq, true, true, true);
    res = HMM.get_training_sequences(seq);
end
