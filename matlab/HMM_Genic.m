% Class for genic HMM sub-model
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko, Jelmer Cnossen, Orr Shomroni
% ------------------------------------------------------------------------
classdef (ConstructOnLoad = true) HMM_Genic < HMM
    properties (GetAccess = protected, SetAccess = protected)
        % Get some statistics, ignore all non-traditional nucleotide bases
        genic_codon_stats
        genic_nucleotide_stats
        genic_nucleotide_count
        genic_nucleotide_probabilities
        intergenic_count_template
        intergenic_count_complement
        coding_codon_count
        intergenic_prob
        center_id
        codon_count
    end
    
    methods (Access = public)
        % Class constructor for Genic HMM sub-model
        % Input:
        %       seq                   - training sequence (with annotated
        %                               genes)
        %       mode                  - 'simple' or 'complex': simple -
        %                               with no nucleotide deletions and
        %                               insertions, complex - with
        %                               nucleotide deletions and insertions
        %       insertion_probability - probability of inserting a
        %                               nucleotide for codon-emitting model
        %       deletion_probability  - probability of deleting a
        %                               nucleotide for codon-emitting model
        %       trainModel            - boolean switch telling whether the
        %                               model should be trained on given
        %                               contig
        % Ouput:
        %       obj                   - HMM object.
        function [obj] = HMM_Genic(seq, mode, insertion_probability, deletion_probability, trainModel)
            if (nargin < 2)
                error('Not enough parameters, consult help.');
            end
            
            if (nargin < 5)
                trainModel = false;
            end
            
            mode = find(strcmp(mode, {'simple', 'complex'}));
            
            if (isempty(mode))
                error('Invalid genic model type');
            end
            
            if (mode == 2 && nargin < 3)
                error('Not enough parameters, consult help.');
            end
            
            seq.Sequence = upper(seq.Sequence);
            
            [obj.genic_codon_stats, obj.genic_nucleotide_stats, obj.GeneCount] = get_statistics(seq, true, true, true, false);
            [obj.intergenic_count_template, obj.intergenic_count_complement] = HMM_Genic.count_intergenic(seq);
            [obj.start_codon_stats, obj.stop_codon_stats] = HMM.get_start_stop_statistics(seq, true, true, true);
            obj.codon_count = HMM_Genic.get_codon_count_in_sequence(seq);
            obj.coding_codon_count = sum(sum(sum(obj.genic_codon_stats)));
            
            obj.genic_nucleotide_count = sum(obj.genic_nucleotide_stats);
            obj.genic_nucleotide_probabilities = obj.genic_nucleotide_stats(1:4) / obj.genic_nucleotide_count;
            
            %fprintf('(%i + %i) / %i\n', obj.intergenic_count_template, obj.intergenic_count_complement, sum(sum(sum(obj.codon_count))));
            obj.intergenic_prob = (obj.intergenic_count_template + obj.intergenic_count_complement) / sum(sum(sum(obj.codon_count)));
            
            % Add center state, which does not emit anything
            obj.add_start;
            obj.center_id = obj.add_state('center', zeros(1, 4));
            obj.add_genic(HMM.Genic_Codons, mode, insertion_probability, deletion_probability);
            obj.Trans(obj.center_id, :) = obj.Trans(obj.center_id, :) / sum(obj.Trans(obj.center_id, :)) * (1 - obj.intergenic_prob);
            obj.add_stop;
            obj.link_stop_start;
            
            obj.fix;
            
            % Now to actually train out algo using Baum-Welsch
            % We have to do it for both template and complementray sequences
            if (trainModel)
                seq = HMM.select_genes(seq, true, true, true);
                obj.train(seq);
            end
        end
    end
    
    methods (Access = private)
        % Builds codon loop and connect it to center
        % Input:
        %       codon                 - char array length 3 (nucleotides)
        %       mode                  - switch (1 for simple model - no
        %                               deletions, 2 for complex model -
        %                               with deletions)
        %       codon_probability     - probability of entering this codon
        %                               loop
        %       insertion_probability - probability of inserting a
        %                               nucleotide for codon-emitting model
        %       deletion_probability  - probability of deleting a
        %                               nucleotide for codon-emitting model
        % Output:
        %       id    - id of newly created codon (begin state)
        function [id] = add_codon_loop(obj, codon, mode, codon_probability, insertion_probability, deletion_probability)
            ids = zeros(1, 5);
            step_probability = 1 - insertion_probability - deletion_probability;
            if (mode == 1)
                step_probability = 1;
            end
            
            % create all nodes and some edges
            ids(1) = obj.add_state([codon '_begin'], zeros(1, 4));
            for aclI = 1:3
                ids(aclI + 1) = obj.add_state([codon, '_', int2str(aclI)], HMM.get_emission_prob_for_nucleotide(codon(aclI)));
                obj.add_edge(ids(aclI), ids(aclI + 1), step_probability);
            end
            ids(5) = obj.add_state([codon '_end'], zeros(1, 4));
            
            if (mode == 1)
                obj.add_edge(ids(4), ids(5), step_probability);
            elseif (mode == 2)
                obj.add_edge(ids(4), ids(5), step_probability + deletion_probability);
            end
            
            % add more edges
            obj.add_edge('center', ids(1), codon_probability); % center -> codon_begin
            obj.add_edge(ids(5), 'center', 1); % codon_end -> center
            if (mode == 2)
                % add edges to account for deletion
                for aclI = 1:3
                    obj.add_edge(ids(aclI), ids(aclI + 2), deletion_probability);
                end
                
                % add edges and nodes to account for insetion before each
                % nucleotide in codon
                for aclI = 1:3
                    id = obj.add_state([codon, '_', int2str(aclI), '_insert'], obj.genic_nucleotide_probabilities);
                    obj.add_edge(ids(aclI), id, insertion_probability);
                    obj.add_edge(id, ids(aclI + 1), 1);
                end
                
                % add edge to account for insertion after last nucleotide
                % in codon
                id = obj.add_state([codon, '_end_insert'], obj.genic_nucleotide_probabilities);
                obj.add_edge(ids(4), id, insertion_probability);
                obj.add_edge(id, ids(5), 1);
            end
            
            id = ids(1);
        end
        
        % Adds genic model into the graph
        % Input:
        %      codons                - a list of genic codons (the ones we
        %                              use to model genic regions;
        %      mode                  - switch (1 for simple model - no
        %                              deletions, 2 for complex model -
        %                              with deletions)
        %      insertion_probability - probability of inserting a
        %                              nucleotide in codon-emitting model
        %      deletion_probability  - probability of deleting a nucleotide
        %                              in codon-emitting model
        function [] = add_genic(obj, codons, mode, insertion_probability, deletion_probability)
            codonsN = length(codons);
            for agI = 1:codonsN
                % Transition probability estimated from the statistics
                p = obj.genic_codon_stats(nt2int(codons{agI}(1)), nt2int(codons{agI}(2)), nt2int(codons{agI}(3))) / obj.coding_codon_count;
                obj.add_codon_loop(codons{agI}, mode, p, insertion_probability, deletion_probability);
            end
        end
        
        % Links stop codons and start codons to center (this is done here
        % because otherwise we can not train our model (hmmtrain wants to
        % start from state number one).
        function [] = link_stop_start(obj)
            % Link start
            obj.add_edge('start_codons_G', 'center', 1);
            % Link stop
%             obj.add_edge('center', 'stop_T', obj.intergenic_prob);
            first = obj.stop_codon_stats(nt2int('T'), nt2int('A'), nt2int('A')) + obj.stop_codon_stats(nt2int('T'), nt2int('G'), nt2int('A'));
            second = obj.stop_codon_stats(nt2int('T'), nt2int('A'), nt2int('G'));
            both = first + second;
            obj.add_edge('center', 'stop_TAA_TGA_1', obj.intergenic_prob * first / both);
            obj.add_edge('center', 'stop_TAG_1', obj.intergenic_prob * second / both);            
        end
    end
    
    methods (Access = private, Static)
        
        % Returns number of codons codon count matrix in a contig (direct
        % and complemnt strand together)
        % Input:
        %      seq - contig
        % Output:
        %      cnt - codon count matrix (as returned by codoncount);
        function [cnt] = get_codon_count_in_sequence(seq)
            [t, s1] = codoncount(seq.Sequence);
            [t, s2] = codoncount(seqrcomplement(seq.Sequence));
            cnt = s1 + s2;
        end
        
        % Returns number of intergenic regions in conting
        % Input:
        %    seq          - contig with annotated genes
        % Output:
        %    directCount  - number of intergenic regions on direcct DNA
        %                   strand.
        %    reverseCount - number of intergenic regions on reverse
        %                   complement DNA strand
        function [directCount, reverseCount] = count_intergenic(seq)
            len = length(seq.gene);
            seqLen = length(seq.Sequence);
            directCount = 0;
            reverseCount = 0;
            
            start = 1;
            for i = 1:len
                if (seq.gene(i).Indices(1) < seq.gene(i).Indices(2))
                    start = seq.gene(i).Indices(2) + 1;
                    directCount = directCount + 1;
                end
            end
            if (start < seqLen)
                directCount = directCount + 1;
            end
            
            start = seqLen;
            for i = len:-1:1
                if (seq.gene(i).Indices(1) > seq.gene(i).Indices(2))
                    start = seq.gene(i).Indices(2) - 1;
                    reverseCount = reverseCount + 1;
                end
            end
            if (start > 1)
                reverseCount = reverseCount + 1;
            end
        end
    end
end