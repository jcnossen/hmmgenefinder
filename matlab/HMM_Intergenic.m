% Base class for intergenic HMMs.
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
classdef (ConstructOnLoad = false) HMM_Intergenic < HMM
    properties (GetAccess = protected, SetAccess = protected)
        intergenic_codon_stats
        intergenic_nucleotide_stats
        intergenic_probabilities
    end
    
    methods (Access = protected)
        % Constructs the whole model. Hope this is a correct one...
        % Input:
        %       prefix - string prefix for all model state names
        %       len    - length of desired intergenic model
        function [] = construct_model(obj, prefix, len)
            top = zeros(1, len);
            middle = zeros(1, len);
            bottom = zeros(1, len);
            beginID = obj.add_state([prefix 'begin'], zeros(1, 4));
            % and begin middle and link it
            for i = 1:len
                top(i) = obj.add_state([prefix 'top_' int2str(i)], zeros(1, 4));
                middle(i) = obj.add_state([prefix 'middle_' int2str(i)], obj.intergenic_probabilities);
                bottom(i) = obj.add_state([prefix 'bottom_' int2str(i)], obj.intergenic_probabilities);
                % bottom -> middle
                obj.add_edge(bottom(i), middle(i), 1 / 3);
                % middle loop
                obj.add_edge(middle(i), middle(i), 1 / 4);
                % middle -> top
                obj.add_edge(middle(i), top(i), 1 / 4);
            end
            
            for i = 1:len - 1
                % middle -> top + 1
                obj.add_edge(middle(i), top(i + 1), 1 / 4);
                % middle -> bottom + 1
                obj.add_edge(middle(i), bottom(i + 1), 1 / 4);
                % bottom -> top + 1
                obj.add_edge(bottom(i), top(i + 1), 1 / 3);
                % top -> bottom + 1
                obj.add_edge(top(i), bottom(i + 1), 1 / 2);
                
                % link top
                obj.add_edge(top(i), top(i + 1), 1 / 2);
                % link bottom
                obj.add_edge(bottom(i), bottom(i + 1), 1 / 3);
            end
            
            endID = obj.add_state([prefix 'end'], zeros(1, 4));
            obj.add_edge(bottom(len), endID, 2 / 3);
            obj.add_edge(middle(len), endID, 1 / 2);
            obj.add_edge(top(len), endID, 1);
            
            middleBeginID = obj.add_state([prefix, 'middle_begin'], obj.intergenic_probabilities);
            obj.add_edge(beginID, middleBeginID, 1 / 3);
            obj.add_edge(beginID, bottom(1), 1 / 3);
            obj.add_edge(beginID, top(1), 1 / 3);
            
            obj.add_edge(middleBeginID, middleBeginID, 1 / 2);
            obj.add_edge(middleBeginID, bottom(1), 1 / 2);
            %obj.add_edge(middleBeginID, top(1), 1 / 3);
        end
    end
       
    methods (Access = public, Static)
        % Selects intergenic regions of appropriate length
        % Input:
        %       seq    - sequence with annotated regions
        %       minLen - minimum length for selected regions
        %       maxLen - maximum length for selected regions
        % Output:
        %       res    - sequence with only those annotations left that
        %                answer length requirements
        function [res] = select_intergenic(seq, minLen, maxLen)
            ind = 0;
            len = length(seq.gene);
            for i = 1:len
                iterLen = abs(seq.gene(i).Indices(1) - seq.gene(i).Indices(2)) + 1;
                if (iterLen >= minLen && iterLen <= maxLen)
                    ind = ind + 1;
                    res.gene(ind) = seq.gene(i);
                end
            end
            res.Sequence = seq.Sequence;
        end
        
        % Constructs intergenic regions from genic regions.
        % Input:
        %       seq               - sequence with annotated genes
        %       withCodons        - boolean switch - should stop and start
        %                           codons be included with the intergenic
        %                           region?
        %       traditionalCodons - boolean switch - select only intergenic
        %                           regions with traditional stop and start
        %                           codons?
        % Output:
        %       res               - sequence with annotated intergenic
        %                           regions, answering given requirements
        function [res] = get_intergenic(seq, withCodons, traditionalCodons)
            len = length(seq.gene);
            seqLen = length(seq.Sequence);
            directCount = 0;
            reverseCount = 0;
            for i = 1:len
                if (seq.gene(i).Indices(1) < seq.gene(i).Indices(2))
                    directCount = directCount + 1;
                    direct(directCount) = seq.gene(i);
                else
                    reverseCount = reverseCount + 1;
                    reverse(reverseCount) = seq.gene(i);
                end
            end
            start = 1;
            ind = 0;
            for i = 1:directCount
                stop = direct(i).Indices(1) - 1;
                if (start <= stop) % No overlapping genes...
                    if (withCodons)
                        if (start > 3 && stop + 3 <= seqLen)
                            stopCodon = seq.Sequence(start - 3:start - 1);
                            startCodon = seq.Sequence(stop + 1:stop + 3);
                            if (~traditionalCodons || (HMM.is_start(startCodon) && HMM.is_stop(stopCodon)))
                                ind = ind + 1;
                                res.gene(ind).Indices(1) = start - 3;
                                res.gene(ind).Indices(2) = stop + 3;
                            end
                        end
                    else
                        ind = ind + 1;
                        res.gene(ind).Indices(1) = start;
                        res.gene(ind).Indices(2) = stop;
                    end
                end
                start = direct(i).Indices(2) + 1;
            end
%             if (start < seqLen)
%                 ind = ind + 1;
%                 res.gene(ind).Indices(1) = start;
%                 res.gene(ind).Indices(2) = seqLen;
%             end
            
            start = seqLen;
            for i = reverseCount:-1:1
                stop = reverse(i).Indices(1) + 1;
                if (start >= stop) % No overlapping genes...
                    if (withCodons)
                        if (stop > 3 && start + 3 <= seqLen)
                            stopCodon = seqrcomplement(seq.Sequence(start + 1:start + 3));
                            startCodon = seqrcomplement(seq.Sequence(stop - 3:stop - 1));
                            if (~traditionalCodons || (HMM.is_start(startCodon) && HMM.is_stop(stopCodon)))
                                ind = ind + 1;
                                res.gene(ind).Indices(1) = start + 3;
                                res.gene(ind).Indices(2) = stop - 3;
                            end
                        end
                    else
                        ind = ind + 1;
                        res.gene(ind).Indices(1) = start;
                        res.gene(ind).Indices(2) = stop;
                    end
                end
                start = reverse(i).Indices(2) - 1;
            end
%             if (start > 1)
%                 ind = ind + 1;
%                 res.gene(ind).Indices(1) = start;
%                 res.gene(ind).Indices(2) = 1;
%             end

            res.Sequence = seq.Sequence;
        end
    end
end