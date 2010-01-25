% Splits given genome in required proportion into training and test
% datasets.
% Input arguments:
%      <f>            - annotated sequence(genome) that requires splitting.
%      [wanted_ratio] - ratio between size of the training set and whole
%                       genome.
%      [impTh]        - Internal parameter - sets number of iterations with
%                       no improvement after which the algorithm
%                       terminates.
% Output arguments:
%      [train]        - Training dataset.
%      [test]         - Testing dataset.
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
function [train, test] = split_contig(f, wanted_ratio, impTh)
    
    % Returns a boolean value - tells us if we can cut the genome at given
    % position.
    function [res] = can_cut(location)
        for pos = 1:n
            gene = f.gene(pos);
            if (min(gene.Indices) <= location && location <= max(gene.Indices))
                res = false;
                return;
            end
        end
        res = true;
    end

    % Returns all genes that start and end in given range
    % g   - genes to select from
    % s,f - start and end of allowed nucleotide range
    function [genes] = get_all_genes(g, s, f)
        genes = [];
        n = length(g.gene);
        for gagI = 1 : n
            if (min(g.gene(gagI).Indices) >= s && max(g.gene(gagI).Indices) <= f)
                genes = [genes g.gene(gagI)];
            end
        end
    end

    % Shifts positions of all genes in f by delta
    function [f] = shift_genes(f, delta)
        n = length(f);
        for sgI = 1:n
            f(sgI).Indices = f(sgI).Indices - delta;
        end
    end

    % Returns number of genes that start and end in given range
    function [cnt] = count_genes(st, fn)
        cnt = 0;
        n = length(f.gene);
        for cgI = 1 : n
            if (min(f.gene(cgI).Indices) >= st && max(f.gene(cgI).Indices) <= fn)
                cnt = cnt + 1;
            end
        end
    end
    
    if (nargin < 1)
        error('Not enough input arguments. Consult help.');
    end
    
    if (nargin < 2)
        wanted_ratio = 0.5;
    end
    
    if (nargin < 3)
        impTh = 5;
    end
    
    seq_length = length(f.Sequence);
    n = length(f.gene);
    
    % For each gene go and try to divide the genome after each gene, make
    % sure we do not cut any genes in the middle and select the ratio
    % closest to 0.5
    best_ratio = Inf;
    best_position = 0;
    
    last_impI = 0;
    last_impJ = 0;
    i = round(n * wanted_ratio) - 1;
    j = round(n * wanted_ratio);
    while ((i > 0 && last_impI < impTh) || (j < n && last_impJ < impTh))
        
        if (i > 0 && last_impI < impTh)
            cur = f.gene(i);
            next = f.gene(i + 1);
            cur_end = max(cur.Indices);
            next_start = min(next.Indices);
            middle = round(mean([cur_end next_start]));
            if (can_cut(middle))
                ratio = count_genes(1, middle) / length(f.gene);
                if (abs(ratio - wanted_ratio) < abs(best_ratio - wanted_ratio))
                    best_ratio = ratio;
                    fprintf('[+] (%i) New best ratio attained - %f\n', i, best_ratio);
                    best_position = middle;
                    last_impI = 0;
                else
                    last_impI = last_impI + 1;
                end
            end
        end
        i = i - 1;
        
        if (j < n && last_impJ < impTh)
            cur = f.gene(j);
            next = f.gene(j + 1);
            cur_end = max(cur.Indices);
            next_start = min(next.Indices);
            middle = round(mean([cur_end next_start]));
            if (can_cut(middle))
                ratio = count_genes(1, middle) / length(f.gene);
                if (abs(ratio - wanted_ratio) < abs(best_ratio - wanted_ratio))
                    best_ratio = ratio;
                    fprintf('[+] (%i) New best ratio attained - %f\n', i, best_ratio);
                    best_position = middle;
                    last_impJ = 0;
                else
                    last_impJ = last_impJ + 1;
                end
            end
        end
        j = j + 1;
    end
    % BTW, this works only coz the genes are sorted in incresing order of
    % their lower index (lower != first)
    
    fprintf('[i] Cutting sequence at %i\n', best_position); 
    train.Sequence = f.Sequence(1:best_position);
    train.gene = get_all_genes(f, 1, best_position);
    test.Sequence = f.Sequence(best_position + 1:seq_length);
    test.gene = get_all_genes(f, best_position + 1, seq_length);
    test.gene = shift_genes(test.gene, best_position);
end