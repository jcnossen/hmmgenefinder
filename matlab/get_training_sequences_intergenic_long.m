% Returns sequences for training HMM (used to export data for GHMM)
% Input:
%       seq - sequence with annotated genes
% Output:
%       res - cell array with intergnenic regions
function [res] = get_training_sequences_intergenic_long(seq)
    % Just in case, but better do it when loading the genome
    seq.Sequence = upper(seq.Sequence);
    
    seq = HMM_Intergenic.get_intergenic(seq, true, true);
    seq = HMM_Intergenic.select_intergenic(seq, 10, Inf);
    seq = HMM.select_genes(seq, true, false, false);
    res = HMM.get_training_sequences(seq);
end
