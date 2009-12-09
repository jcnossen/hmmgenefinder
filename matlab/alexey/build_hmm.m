function [trans, emit, symbols, stateNames] = build_hmm(seq, insertion_probability, deletion_probability)
    % Some global variables
    stateCount = 0;
    trans = [];
    emit = zeros(0, 4);
    stateNames = cell(1, 0);
    
    % Predefine some stuff
    genic_codon_stats = [];
    genic_codon_count = 0;
    genic_nucleotide_stats = [];
    genic_nucleotide_count = 0;
    genic_nucleotide_probabilities = [];
    genic_codons = {'AAA', 'AAC', 'AAT', 'AAG', 'ACA', 'ACC', 'ACT', 'ACG', 'ATA', 'ATC', 'ATT', 'AGA', 'AGC', 'AGT', 'AGG', 'CAA', 'CAC', 'CAT', 'CAG', 'CCA', 'CCC', 'CCT', 'CCG', 'CTA', 'CTC', 'CTT', 'CTG', 'CGA', 'CGC', 'CGT', 'CGG', 'TAC', 'TAT', 'TCA', 'TCC', 'TCT', 'TCG', 'TTA', 'TTC', 'TTT', 'TTG', 'TGC', 'TGT', 'TGG', 'GAA', 'GAC', 'GAT', 'GAG', 'GCA', 'GCC', 'GCT', 'GCG', 'GTA', 'GTC', 'GTT', 'GGA', 'GGC', 'GGT', 'GGG'};
    
    % Making life easy
    if (nargin < 1)
        error('Not enough parameters, consult help.');
    end
    
    if (nargin < 2)
        insertion_probability = 10 ^ -8;
        deletion_probability = 10 ^ -8;
        fprintf('[!] No insertion and deletion probabilities defined.\n');
        fprintf('[i] Assuming Pins = %e, Pdel = %e\n', insertion_probability, deletion_probability);
    end
    
    %----------------------------------------------------------------------
    % Some support functions to make life easier
    
    % Return emition probabilities for nucleotide
    % nuc  - nucleotide given as number
    % prob - emition probabilities
    function [p] = get_emition_prob_for_nucleotide(nuc)
        p = zeros(1, 4);
        nucNum = nt2int(nuc);
        p(nucNum) = 1;
    end
    
    % Adds a new state with given name
    % name - state name
    % row  - id of newly added state
    % prob - emition probability
    function [row] = add_state(name, prob)
        if (nargin < 2)
            prob = zeros(1, 4);
        end
        stateCount = stateCount + 1;
        trans = [trans, zeros(stateCount - 1, 1); zeros(1, stateCount)];
        emit = [emit; prob];
        stateNames{stateCount} = name;
        
        row = stateCount;
    end

    % Adds a new edge between two sates
    % i,j - state ids or state names
    % p   - probability of transition
    function [] = add_edge(i, j, p)
        if (isa(i, 'char'))
            i = get_id_by_state_name(i);
        end
        if (isa(j, 'char'))
            j = get_id_by_state_name(j);
        end
        trans(i, j) = p;
    end

    % Build codon loop and connect it to center
    % codon - char array length 3 (nucleotides)
    % id    - id of newly created codon (begin state)
    function [id] = add_codon_loop(codon, codon_probability, insertion_probability, deletion_probability)
        ids = zeros(1, 5);
        step_probability = 1 - insertion_probability - deletion_probability;
        
        % create all nodes and some edges
        ids(1) = add_state([codon '_begin'], [0 0 0 0]);
        for aclI = 1:3
            ids(aclI + 1) = add_state([codon, '_', codon(aclI)], get_emition_prob_for_nucleotide(codon(aclI)));
            add_edge(ids(aclI), ids(aclI + 1), step_probability);
        end
        ids(5) = add_state([codon '_end'], [0 0 0 0]);
        
        % add more edges
        add_edge('center', ids(1), codon_probability); % center -> codon_begin
        add_edge(ids(5), 'center', 1); % codon_end -> center
        % add edges to account for deletion
        for aclI = 1:3
            add_edge(ids(aclI), ids(aclI + 2), deletion_probability);
        end
        % add edges and nodes to account for insetion before each
        % nucleotide in codon
        for aclI = 1:3
            id = add_state([codon, '_', codon(aclI), '_insert'], genic_nucleotide_probabilities);
            add_edge(ids(aclI), id, insertion_probability);
            add_edge(id, ids(aclI + 1), 1);
        end
        % add edge to account for insertion after last nucleotide in codon
        id = add_state([codon, '_end_insert'], genic_nucleotide_probabilities);
        add_edge(ids(4), id, insertion_probability);
        add_edge(id, ids(5), 1);
        
        id = ids(1);
    end

    % Returns state id from a state name
    function [id] = get_id_by_state_name(name)
        id = 0;
        for idI = 1:stateCount
            if (strcmpi(stateNames{idI}, name))
                id = idI;
                return;
            end
        end
    end

    % Adds intergenic model into global model
    function [] = add_genic(codons, insertion_probability, deletion_probability)
        genic_nucleotide_probabilities = genic_nucleotide_stats(1:4) / genic_nucleotide_count;
        codonsN = length(codons);
        for agI = 1:codonsN
            % Transition probability estimated from the statistics
            p = genic_codon_stats(nt2int(codons{agI}(1)), nt2int(codons{agI}(2)), nt2int(codons{agI}(3))) / genic_codon_count;
            add_codon_loop(codons{agI}, p, insertion_probability, deletion_probability);
        end
    end
    %----------------------------------------------------------------------

    % Get some statistics, ignore all non-traditional nucleotide bases
    [genic_codon_stats, genic_nucleotide_stats] = get_genic_statistics(seq, true);
    genic_codon_count = sum(sum(sum(genic_codon_stats)));
    genic_nucleotide_count = sum(genic_nucleotide_stats);
    
    % Add center state, which does not emit anything
    add_state('center', [0 0 0 0]);
    add_genic(genic_codons, insertion_probability, deletion_probability);
    %add_stop_codons;
    %add_intergenic;
    %add_start_codons;
    
    % Now to actually train out algo using Baum-Welsch
    % We have to do it for both template and complementray sequences
    %hmmtrain(sequence, trans, emit, 'Algorithm', 'BaumWelch', 'Verbose', true);
    
    symbols = {'A', 'C', 'T', 'G'};
end