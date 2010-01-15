% HMM class for simple intergenic model.
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko, Jelmer Cnossen, Orr Shomroni
% ------------------------------------------------------------------------
classdef (ConstructOnLoad = true) HMM_Intergenic_Dumb < HMM_Intergenic
    methods (Access = public)
        % Class constructor.
        % Input:
        %       seq        - sequence with annotaged genes
        %       trainModel - boolean switch telling whether the model
        %                    should be trained
        % Output:
        %       obj        - HMM object
        function [obj] = HMM_Intergenic_Dumb(seq, trainModel)
            if (nargin < 2)
                trainModel = false;
            end
            seq.Sequence = upper(seq.Sequence);
            [obj.start_codon_stats, obj.stop_codon_stats] = HMM.get_start_stop_statistics(seq, true, true, true);
            
            seq = HMM_Intergenic.get_intergenic(seq, true, true);
            avg_len = HMM_Intergenic_Dumb.get_average_length(seq);
            [obj.intergenic_codon_stats, obj.intergenic_nucleotide_stats, obj.GeneCount] = get_statistics(seq, true, false, false, false);
            obj.intergenic_probabilities = obj.intergenic_nucleotide_stats(1:4) / sum(obj.intergenic_nucleotide_stats(1:4));
            
            obj.add_stop;
            id = obj.add_state('intergenic_dumb', obj.intergenic_probabilities);
            obj.add_edge(id, id, 1 - 1 / avg_len);
            obj.add_start;
            
            % Link stop and start to intergenic model
            obj.add_edge('stop_TAA_TGA_3', 'intergenic_dumb', 1);
            obj.add_edge('stop_TAG_3', 'intergenic_dumb', 1);
            obj.add_edge('intergenic_dumb', 'start_codons_AGT', 1 / avg_len);
            
            if (trainModel)
                obj.train(seq);
            end
        end
    end
    
    methods (Access = private, Static)
        % Calculated average length of intergenic sequences
        % Input:
        %       seq - sequence with annotated regions
        % Output:
        %       len - average length of intergenic seuqnces
        function [len] = get_average_length(seq)
            len = 0;
            cnt = length(seq.gene);
            for i = 1:cnt
                len = len + abs(seq.gene(i).Indices(1) - seq.gene(i).Indices(2)) + 1;
            end
            len = floor(len / cnt);
        end
    end
end