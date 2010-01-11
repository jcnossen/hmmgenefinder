% HMM class for short complex intergenic model.
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko, Jelmer Cnossen, Orr Shomroni
% ------------------------------------------------------------------------
classdef (ConstructOnLoad = true) HMM_Intergenic_Short < HMM_Intergenic
    methods (Access = public)
        % Class constructor.
        % Input:
        %       seq        - sequence with annotaged genes
        %       trainModel - boolean switch telling whether the model
        %                    should be trained
        % Output:
        %       obj        - HMM object
        function [obj] = HMM_Intergenic_Short(seq, trainModel)
            if (nargin < 2)
                trainModel = false;
            end
            seq.Sequence = upper(seq.Sequence);
            
            seq = HMM_Intergenic.get_intergenic(seq, true, true);
            seq = HMM_Intergenic.select_intergenic(seq, 1, 14);
            [obj.intergenic_codon_stats, obj.intergenic_nucleotide_stats, obj.GeneCount] = get_statistics(seq, true, false, false, true);
            obj.intergenic_probabilities = obj.intergenic_nucleotide_stats(1:4) / sum(obj.intergenic_nucleotide_stats(1:4));
            
            obj.construct_model('intergenic_short_', 9);
            if (trainModel)
                obj.train(seq);
            end
        end
    end
end