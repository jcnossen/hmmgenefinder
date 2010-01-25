% HMM class for short complex intergenic model.
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
classdef (ConstructOnLoad = true) HMM_Intergenic_Short < HMM_Intergenic
    methods (Access = public)
        % Class constructor.
        % Input:
        %       seq        - sequence with annotated genes
        %       trainModel - boolean switch telling whether the model
        %                    should be trained
        % Output:
        %       obj        - HMM object
        function [obj] = HMM_Intergenic_Short(seq, trainModel)
            if (nargin < 2)
                trainModel = false;
            end
            seq.Sequence = upper(seq.Sequence);
            [obj.start_codon_stats, obj.stop_codon_stats] = HMM.get_start_stop_statistics(seq, true, true, true);
            
            seq = HMM_Intergenic.get_intergenic(seq, true, true);
            seq = HMM_Intergenic.select_intergenic(seq, 1, 14);
            [obj.intergenic_codon_stats, obj.intergenic_nucleotide_stats, obj.GeneCount] = get_statistics(seq, true, false, false, false);
            obj.intergenic_probabilities = obj.intergenic_nucleotide_stats(1:4) / sum(obj.intergenic_nucleotide_stats(1:4));
            
            obj.add_stop;
            obj.construct_model('intergenic_short_', 9);
            obj.add_start;
            
            % Link stop and start to intergenic model
            obj.add_edge('stop_TAA_TGA_3', 'intergenic_short_begin', 1);
            obj.add_edge('stop_TAG_3', 'intergenic_short_begin', 1);
            obj.add_edge('intergenic_short_end', 'start_codons_AGT', 1);
            
            obj.fix;
            
            if (trainModel)
                seq = HMM.select_genes(seq, true, false, false);
                obj.insert_state(1, 'entry_point', zeros(1, 4));
                obj.add_edge(1, 'stop_T', 1);
                obj.train(seq);
                obj.delete_state('entry_point');
            end
        end
    end
end