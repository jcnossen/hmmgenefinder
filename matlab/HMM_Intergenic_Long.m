% HMM class for long complex intergenic model.
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
classdef (ConstructOnLoad = true) HMM_Intergenic_Long < HMM_Intergenic
    methods (Access = public)
        % Class constructor.
        % Input:
        %       seq        - sequence with annotaged genes
        %       trainModel - boolean switch telling whether the model
        %                    should be trained
        % Output:
        %       obj        - HMM object
        function [obj] = HMM_Intergenic_Long(seq, trainModel)
            if (nargin < 2)
                trainModel = false;
            end
            seq.Sequence = upper(seq.Sequence);
            
            seq = HMM_Intergenic.get_intergenic(seq, true, true);
            seq = HMM_Intergenic.select_intergenic(seq, 10, Inf);
            [obj.intergenic_codon_stats, obj.intergenic_nucleotide_stats, obj.GeneCount] = get_statistics(seq, true, false, false, true);
            obj.intergenic_probabilities = obj.intergenic_nucleotide_stats(1:4) / sum(obj.intergenic_nucleotide_stats(1:4));
            
            obj.construct_model('intergenic_long_', 44);
            if (trainModel)
                obj.train(seq);
            end
        end
    end
end