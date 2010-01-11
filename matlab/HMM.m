% Base class for HMM models
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko, Jelmer Cnossen, Orr Shomroni
% ------------------------------------------------------------------------
classdef (ConstructOnLoad = false) HMM < handle
    properties
        StateCount = 0;
        GeneCount = 0;
        Trans = [];
        Emit = zeros(0, 4);
        StateNames = cell(1, 0);
        Symbols = num2cell(int2nt(1:4));
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        start_codon_stats
        stop_codon_stats
    end

    properties (SetAccess = private, Constant)
        % Predefine some stuff
        Genic_Codons = {'AAA', 'AAC', 'AAT', 'AAG', 'ACA', 'ACC', 'ACT', 'ACG', 'ATA', 'ATC', 'ATT', 'ATG', 'AGA', 'AGC', 'AGT', 'AGG', 'CAA', 'CAC', 'CAT', 'CAG', 'CCA', 'CCC', 'CCT', 'CCG', 'CTA', 'CTC', 'CTT', 'CTG', 'CGA', 'CGC', 'CGT', 'CGG', 'TAC', 'TAT', 'TCA', 'TCC', 'TCT', 'TCG', 'TTA', 'TTC', 'TTG', 'TTT', 'TGC', 'TGT', 'TGG', 'GAA', 'GAC', 'GAT', 'GAG', 'GCA', 'GCC', 'GCT', 'GCG', 'GTA', 'GTC', 'GTT', 'GTG', 'GGA', 'GGC', 'GGT', 'GGG'};
        Start_Codons = {'ATG', 'GTG', 'TTG'};
        Stop_Codons  = {'TAA', 'TGA', 'TAG'};
    end
    
    methods (Access = public)
        % HMM object copy constructor.
        % Creates an HMM model identical to the one given as a parameter.
        % Input:
        %       model - HMM object to create a copy of.
        % Output:
        %       obj   - HMM object copy.
        function [obj] = HMM(model)
            if (nargin > 0)
                obj.StateCount = model.StateCount;
                obj.GeneCount = model.GeneCount;
                obj.Trans = model.Trans;
                obj.Emit = model.Emit;
                obj.StateNames = model.StateNames;
                obj.Symbols = model.Symbols;
            end
        end
    end
    
    methods (Access = public, Static)
        % Returns true if the codon is a (traditional) start codon
        % Input:
        %       codon - a string of 3 nucleotide letters
        % Output:
        %       res   - boolean value with result.
        function [res] = is_start(codon)
            found = find(strcmpi(codon, HMM.Start_Codons), 1);
            res = ~isempty(found);
        end
        
        % Returns true is the codon is a (traditions) stop codon
        % Input:
        %       codon - a string of 3 nucleotide letters
        % Output:
        %       res   - boolean value with result.
        function [res] = is_stop(codon)
            found = find(strcmpi(codon, HMM.Stop_Codons), 1);
            res = ~isempty(found);
        end
    end
    
    methods
        % Adds a new state with given name
        % Input:
        %       name - name of the new state
        %       prob - emission probabilities for new state (array of
        %              length 4)
        % Output:
        %       row  - id of newly added state
        function [row] = add_state(obj, name, prob)
            if (nargin < 2)
                prob = zeros(1, 4);
            end
            obj.StateCount = obj.StateCount + 1;
            obj.Trans = [obj.Trans, zeros(obj.StateCount - 1, 1); zeros(1, obj.StateCount)];
            obj.Emit = [obj.Emit; prob];
            obj.StateNames{obj.StateCount} = name;
            
            row = obj.StateCount;
        end
        
        % Inserts a new sate in given position
        % Input:
        %      pos  - position for insertion
        %      name - name of state
        %      prob - emission probabilities for new state (array of length
        %             4)
        % Output:
        %      res  - boolean value representing success of operation
        function [res] = insert_state(obj, pos, name, prob)
            res = false;
            if (pos > 0 && pos <= obj.StateCount)
                for i = obj.StateCount:-1:pos
                    obj.StateNames{i + 1} = obj.StateNames{i};
                end
                obj.StateNames{pos} = name;
                obj.Emit = [obj.Emit(1:pos - 1, :); prob; obj.Emit(pos:obj.StateCount, :)];
                obj.Trans = [obj.Trans(1:pos - 1, 1:pos - 1), zeros(pos - 1, 1), obj.Trans(1:pos - 1, pos:obj.StateCount); zeros(1, obj.StateCount + 1); obj.Trans(pos:obj.StateCount, 1:pos - 1), zeros(obj.StateCount - pos + 1, 1), obj.Trans(pos:obj.StateCount, pos:obj.StateCount)];
                obj.StateCount = obj.StateCount + 1;
                res = true;
            end
        end
        
        % Adds a new edge between two sates
        % Input:
        %       i, j - state ids or state names to define transition i -> j
        %       p    - probability of the transition
        function [] = add_edge(obj, i, j, p)
            %fprintf('%s -> %s\n', i, j);
            if (isa(i, 'char'))
                i = obj.get_id_by_state_name(i);
            end
            if (isa(j, 'char'))
                j = obj.get_id_by_state_name(j);
            end
            obj.Trans(i, j) = p;
        end
        
        % Trains the model
        % Input:
        %       seq - training sequenced (with annotated regions)
        function [] = train(obj, seq, varargin)
            pemit = zeros(obj.StateCount, length(obj.Symbols));
            ptrans = zeros(obj.StateCount);
            algo = 'BaumWelch';
            if (nargin > 2)
                for j = 1:2:nargin - 2
                    pname = varargin{j};
                    pval = varargin{j + 1};
                    k = strmatch(lower(pname), {'pseudoemissions', 'pseudotransitions'});
                    if (~isempty(k))
                        if (k == 1)
                            pemit = pval;
                        elseif (k == 2)
                            ptrans = pval;
                        end
                        algo = 'Viterbi';
                    end
                end
            end
            
            seq = HMM.get_training_sequences(seq);
            fprintf('[i] Training model on %i sequences.\n', length(seq));
            [obj.Trans, obj.Emit] = hmmtrain(seq, obj.Trans, obj.Emit, 'Algorithm', algo, 'Verbose', true, 'Symbols', obj.Symbols, 'Pseudoemissions', pemit, 'Pseudotransitions', ptrans);
        end
        
        % Returns most probable path for given sequence
        % Input:
        %       seq         - a string of nucleoties that is fed to Viterbi
        %                     algorithm
        % Output:
        %       path        - array of states constituting the most
        %                     probable path for the sequence
        function [path] = viterbi(obj, seq)
            seq = num2cell(seq);
            path = my_hmmviterbi(seq, obj.Trans, obj.Emit, 'Statenames', obj.StateNames, 'Symbols', obj.Symbols);
        end
        
        % Removes fields from HMM that do no emit anything. This is
        % required for Viterbi.
        % Also normalizes transition and emission probabilities to sum up
        % to 1.
        function [] = fix(obj)
            noEmitArr = [];
            for i = 1:obj.StateCount
                if (sum(obj.Emit(i, :)) == 0)
                    noEmitArr = [noEmitArr i];
                end
            end
            
            n = length(noEmitArr);
            for i = 1:n
                noEmit = noEmitArr(i) - i + 1;
                fprintf('[i] Removing node %s.\n', obj.StateNames{noEmit});
                from = [];
                to = [];
                for j = 1:obj.StateCount
                    if (obj.Trans(j, noEmit) > 0)
                        from = [from j];
                    end
                    if (obj.Trans(noEmit, j) > 0)
                        to = [to j];
                    end
                end
                
                fromN = length(from);
                toN = length(to);
                for j = 1:fromN
                    for k = 1:toN
                        obj.add_edge(from(j), to(k), obj.Trans(from(j), to(k)) + obj.Trans(from(j), noEmit) * obj.Trans(noEmit, to(k)));
                    end
                end
                
                allowed = setdiff(1:obj.StateCount, noEmit);
                obj.Trans = obj.Trans(allowed, allowed);
                obj.Emit = obj.Emit(allowed, :);
                obj.StateCount = obj.StateCount - 1;
                obj.StateNames = obj.StateNames(allowed);
            end
            
            for i = 1:obj.StateCount
                s = sum(obj.Trans(i, :));
                if (s ~= 1 && s > 0)
                    obj.Trans(i, :) = obj.Trans(i, :) / s;
                end
                
                s = sum(obj.Emit(i, :));
                if (s ~= 1 && s > 0)
                    obj.Emit(i, :) = obj.Emit(i, :) / s;
                end
            end
        end
    end
    
    methods (Access = protected)
        % Returns state id from a state name
        % Input:
        %       name - state name
        % Output:
        %       id   - corresponding state ID.
        function [id] = get_id_by_state_name(obj, name)
            id = 0;
            for idI = 1:obj.StateCount
                if (strcmpi(obj.StateNames{idI}, name))
                    id = idI;
                    return;
                end
            end
        end
        
        % Adds stop codons to the model
        function [] = add_stop(obj)
            id1(1) = obj.add_state('stop_TAA_TGA_1', HMM.get_emission_prob_for_nucleotide('T'));
            count_TAA = obj.stop_codon_stats(nt2int('T'), nt2int('A'), nt2int('A'));
            count_TGA = obj.stop_codon_stats(nt2int('T'), nt2int('G'), nt2int('A'));
            count_both = count_TAA + count_TGA;
            id1(2) = obj.add_state('stop_TAA_TGA_2', (HMM.get_emission_prob_for_nucleotide('A') * count_TAA + HMM.get_emission_prob_for_nucleotide('G') * count_TGA) ./ count_both);
            id1(3) = obj.add_state('stop_TAA_TGA_3', HMM.get_emission_prob_for_nucleotide('A'));
            obj.add_edge(id1(1), id1(2), 1);
            obj.add_edge(id1(2), id1(3), 1);
            
            id2(1) = obj.add_state('stop_TAG_1', HMM.get_emission_prob_for_nucleotide('T'));
            id2(2) = obj.add_state('stop_TAG_2', HMM.get_emission_prob_for_nucleotide('A'));
            id2(3) = obj.add_state('stop_TAG_3', HMM.get_emission_prob_for_nucleotide('G'));
            
            obj.add_edge(id2(1), id2(2), 1);
            obj.add_edge(id2(2), id2(3), 1);

%              % Single entry point to stop codons!
%              id_main = obj.add_state('stop_T', HMM.get_emission_prob_for_nucleotide('T'));
%               
%              count_TAA = obj.stop_codon_stats(nt2int('T'), nt2int('A'), nt2int('A'));
%              count_TGA = obj.stop_codon_stats(nt2int('T'), nt2int('G'), nt2int('A'));
%              count_both = count_TAA + count_TGA;
%              id1(1) = obj.add_state('stop_TAA_TGA_2', (HMM.get_emission_prob_for_nucleotide('A') * count_TAA + HMM.get_emission_prob_for_nucleotide('G') * count_TGA) ./ count_both);
%              id1(2) = obj.add_state('stop_TAA_TGA_3', HMM.get_emission_prob_for_nucleotide('A'));
%              obj.add_edge(id1(1), id1(2), 1);
%              
%              count_TAG = obj.stop_codon_stats(nt2int('T'), nt2int('A'), nt2int('G'));
%              id2(1) = obj.add_state('stop_TAG_2', HMM.get_emission_prob_for_nucleotide('A'));
%              id2(2) = obj.add_state('stop_TAG_3', HMM.get_emission_prob_for_nucleotide('G'));
%              obj.add_edge(id2(1), id2(2), 1);
%               
%              % Now to link single entry point to 2 stop codon models
%              count_all = count_both + count_TAG;
%              obj.add_edge(id_main, id1(1), count_both / count_all);
%              obj.add_edge(id_main, id2(1), count_TAG / count_all);
        end
        
        % Adds start codons to the model
        function [] = add_start(obj)
            count_ATG = obj.start_codon_stats(nt2int('A'), nt2int('T'), nt2int('G'));
            count_GTG = obj.start_codon_stats(nt2int('G'), nt2int('T'), nt2int('G'));
            count_TTG = obj.start_codon_stats(nt2int('T'), nt2int('T'), nt2int('G'));
            count_all = count_ATG + count_GTG + count_TTG;
            id(1) = obj.add_state('start_codons_AGT', (HMM.get_emission_prob_for_nucleotide('A') * count_ATG + HMM.get_emission_prob_for_nucleotide('T') * count_TTG + HMM.get_emission_prob_for_nucleotide('G') * count_GTG) ./ count_all);
            id(2) = obj.add_state('start_codons_T', HMM.get_emission_prob_for_nucleotide('T'));
            id(3) = obj.add_state('start_codons_G', HMM.get_emission_prob_for_nucleotide('G'));
            obj.add_edge(id(1), id(2), 1);
            obj.add_edge(id(2), id(3), 1);
        end
    end
    
    methods (Access = protected, Static)
        
        % Returns only those genes, that fulfil the requirements.
        % Input:
        %       seq                    - sequence with annotated regions
        %       traditionalNucleotides - boolean switch to return only
        %                                genes with traditional nucleotides
        %       traditionalCodons      - boolean switch to return only
        %                                genes with traditional start and
        %                                stop codons
        %       length3                - boolean switch to return only
        %                                genes with length devisible by 3
        % Output:
        %       res                    - same sequence with annotations
        %                                kept only for those genes that
        %                                fulfill the requirements
        function [res] = select_genes(seq, traditionalNucleotides, traditionalCodons, length3)
            n = length(seq.gene);
            geneC = 0;
            for i = 1:n
                s = min(seq.gene(i).Indices);
                f = max(seq.gene(i).Indices);
                len = f - s + 1;
                if (~length3 || mod(len, 3) == 0)
                    str = seq.Sequence(s:f);
                    if (seq.gene(i).Indices(1) > seq.gene(i).Indices(2))
                        str = seqrcomplement(str);
                    end
                    if (len >= 3)
                        startCodon = str(1:3);
                        stopCodon = str(len - 2:len);
                    else
                        startCodon = '';
                        stopCodon = '';
                    end
                    if (~traditionalCodons || (HMM.is_start(startCodon) && HMM.is_stop(stopCodon)))
                        int_str = nt2int(str);
                        if (~traditionalNucleotides || sum(int_str > 4) == 0)
                            geneC = geneC + 1;
                            res.gene(geneC).Indices = seq.gene(i).Indices;
                        end
                    end
                end
            end
            res.Sequence = seq.Sequence;
        end
        
        % Return codon emission count for codons in cell array.
        % Input:
        %       stats  - array with codon statistics (as returned by
        %                codoncount function)
        %       codons - a cell array of codons count should be returned
        %                for
        % Output:
        %       p      - array of counts for requested codons.
        function [p] = get_codon_count(stats, codons)
            n = length(codons);
            p = zeros(1, n);
            for i = 1:n
                codon = codons{i};
                p(i) = stats(nt2int(codon(1)), nt2int(codon(2)), nt2int(codon(3)));
            end
        end
        
        % Return emition probabilities for nucleotide
        % Input:
        %       nuc  - nucleotide given as number (nt2int function).
        % Output:
        %       p    - an array with emission probabilities, that emits
        %              given codon with probability 1
        function [p] = get_emission_prob_for_nucleotide(nuc)
            p = zeros(1, 4);
            nucNum = nt2int(nuc);
            p(nucNum) = 1;
        end
        
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

        % Returns start codon statistics for given sequence
        % Input:
        %       seq                    - sequence with annotated genes
        %       traditionalNucleotides - boolean switch to return only
        %                                genes with traditional nucleotides
        %       traditionalCodons      - boolean switch to return only
        %                                genes with traditional start and
        %                                stop codons
        %       length3                - boolean switch to return only
        %                                genes with length devisible by 3
        % Output:
        %       startC                 - array with codon statistics for
        %                                start codons
        %       stopC                  - array with codon statistics for
        %                                stop codons
        function [startC, stopC] = get_start_stop_statistics(seq, traditionalNucleotides, length3, traditionalCodons)
            startC = zeros([16 16 16]);
            stopC = zeros([16 16 16]);
            n = length(seq.gene);
            for i = 1:n
                st = seq.gene(i).Indices(1);
                fn = seq.gene(i).Indices(2);
                len = abs(fn - st) + 1;
                if (~length3 || mod(len, 3) == 0)
                    rev = false;
                    if (st < fn)
                        start = seq.Sequence(st:st + 2);
                        stop = seq.Sequence(fn - 2:fn);
                        gene = seq.Sequence(st:fn);
                    else
                        start = seq.Sequence(st - 2:st);
                        stop = seq.Sequence(fn:fn + 2);
                        gene = seqrcomplement(seq.Sequence(fn:st));
                        rev = true;
                    end
                    
                    if (~traditionalNucleotides || sum(nt2int(gene) > 4) == 0)
                        if (rev == true)
                            start = seqrcomplement(start);
                            stop = seqrcomplement(stop);
                        end
                        if (~traditionalCodons || (HMM.is_start(start) && HMM.is_stop(stop)))
                            startC(nt2int(start(1)), nt2int(start(2)), nt2int(start(3))) = startC(nt2int(start(1)), nt2int(start(2)), nt2int(start(3))) + 1;
                            stopC(nt2int(stop(1)), nt2int(stop(2)), nt2int(stop(3))) = stopC(nt2int(stop(1)), nt2int(stop(2)), nt2int(stop(3))) + 1;
                        end
                    end
                end
            end
        end
        
    end
end