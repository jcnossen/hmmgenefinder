function [train, test] = split_ecoli(name, wanted_ratio, impTh)
    
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

    % Returns all genes who start and end in given range
    % g   - genes to select from
    % s,f - start and end of allowed nucleotide range
    function [genes] = get_all_genes(g, s, f)
        genes = [];
        n = length(g.gene);
        for gagI = 1 : n
            if (min(g.gene(gagI).Indices) > s && max(g.gene(gagI).Indices) < f)
                genes = [genes g.gene(gagI)];
            end
        end
    end
    
    if (nargin < 1)
        name = 'AE005174';
    end
    
    if (nargin < 2)
        wanted_ratio = 1;
    end
    
    if (nargin < 2)
        impTh = 5;
    end
    
    g = get_bank_file(name);
    f = featuresparse(g.Features);
    
    seq_length = length(g.Sequence);
    n = length(f.gene);
    
    % For each gene go and try to divide the genome after each gene, make
    % sure we do not cut any genes in the middle and select the ratio
    % closest to 0.5
    best_ratio = Inf;
    best_position = 0;
    
    last_impI = 0;
    last_impJ = 0;
    i = round(n / 2) - 1;
    j = round(n / 2);
    while ((i > 0 && last_impI < impTh) || (j < n && last_impJ < impTh))
        
        if (i > 0 && last_impI < impTh)
            cur = f.gene(i);
            next = f.gene(i + 1);
            cur_end = max(cur.Indices);
            next_start = min(next.Indices);
            middle = round(mean([cur_end next_start]));
            if (can_cut(middle))
                ratio = middle / (seq_length - middle);
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
                ratio = middle / (seq_length - middle);
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
    
    train.Sequence = g.Sequence(1:best_position);
    train.gene = get_all_genes(f, 1, best_position);
    test.Sequence = g.Sequence(best_position + 1:seq_length);
    test.gene = get_all_genes(f, best_position + 1, seq_length);
end