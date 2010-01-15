% Performs perfomance assesment of given HMM on given testing contig.
% It runs Viterbi algorithm over the whole contig sequence, finds genes in
% results returned by Viterbi and then tries to match genes it got with
% genes annotated in the contig.
% Input:
%       <hmm>     - HMM object representing HMM we would like to test.
%       <seq>     - sequence with annotated genes.
% Output:
%       <testRes> - structure with test results
function [testRes] = test_hmm(hmm, seq)

    % Finds genes in states array returned by Viterbi
    % Input:
    %       <genes>   - array of annotated gene we have to keep filling
    %       <states>  - cell array of states as returned by Viterbi
    %       <reverse> - boolean switch telling whether the genes are on
    %                   complementary contig (used to swap indices when
    %                   adding to array)
    % Output:
    %       <genes>   - array with annotated genes
    function [genes] = get_genes(genes, states, reverse)
        n = length(states);
        start = 0;
        finish = 0;
        geneC = length(genes);
        for iGG = 1:n
            if (~strncmp('intergenic', states{iGG}, 10))
                % Got a gene!
                if (start == 0)
                    start = iGG;
                else
                    finish = iGG;
                end
            elseif (start ~= 0)
                geneC = geneC + 1;
                if (reverse)
                    genes(geneC).Indices(1) = finish;
                    genes(geneC).Indices(2) = start;
                else
                    genes(geneC).Indices(1) = start;
                    genes(geneC).Indices(2) = finish;
                end
                start = 0;
            end
        end
    end
    
    % Replaces all non-traditional nucleotides in given sequence by
    % randomly choosing one of the possibilities.
    % Input:
    %       <seq> - string representing DNA sequence
    % Output:
    %       <seq> - same sequence but now without non-traditional
    %               nucleotides.
    function [seq] = fix_sequence(seq)
        n = length(seq);
        
        seq = upper(seq);
        for i = 1:n
            if (nt2int(seq(i)) > 4)
                switch (seq(i))
                    case {'r', 'R'}
                        src = ['G', 'A'];
                    case {'y', 'Y'}
                        src = ['T', 'C'];
                    case {'k', 'K'}
                        src = ['G', 'T'];
                    case {'m', 'M'}
                        src = ['A', 'C'];
                    case {'s', 'S'}
                        src = ['G', 'C'];
                    case {'w', 'W'}
                        src = ['A', 'T'];
                    case {'b', 'B'}
                        src = ['G', 'T', 'C'];
                    case {'d', 'D'}
                        src = ['G', 'A', 'T'];
                    case {'h', 'H'}
                        src = ['A', 'C', 'T'];
                    case {'v', 'V'}
                        src = ['G', 'C', 'A'];
                    case {'n', 'N'}
                        src = ['A', 'G', 'C', 'T'];
                    case {'-', ' ', '*'}
                        error('Gaps are not allowed in sequences');
                end
                seq(i) = char(randsrc(1, 1, src + 0));
            end
        end
    end

    % Takes two segments given as 2-element arrays and returns their
    % intersection as a segment.
    % Input:
    %       <a, b>  - segments to return intersection of
    % Output:
    %        <l, r> - left and right boundaries of intersection segment
    function [l, r] = get_intersection(a, b)
        if (a(1) > a(2))
            t = a(1);
            a(1) = a(2);
            a(2) = t;
        end
        if (b(1) > b(2))
            t = b(1);
            b(1) = b(2);
            b(2) = t;
        end
        
        l = max([a(1), b(1)]);
        r = min([a(2), b(2)]);
    end

    % An algorithm for matching genes we got with those annotated in the
    % contig. Basically we have to solve assignment problem here.
    % Input:
    %       <a>   - Genes annotated in the contig
    %       <b>   - Genes we got from Viterbi
    % Output:
    %       <res> - Structure containing results of our matching with
    %               scores and counts.
    function [res] = match_genes(a, b)
        m = length(a);
        n = length(b);
        
        srows = zeros(m * n, 3);
        rowsC = 0;
        usedA = false(1, m);
        usedB = false(1, n);
        
        for i = 1:m
            % each gene in a
            for j = 1:n
                minS = min([a(i).Indices, b(j).Indices]);
                maxS = max([a(i).Indices, b(j).Indices]);
                len = maxS - minS + 1;
                [l, r] = get_intersection(a(i).Indices, b(j).Indices);
                overlap = r - l + 1;
                
                if (overlap > 0)
                    rowsC = rowsC + 1;
                    srows(rowsC, 1) = overlap / len;
                    srows(rowsC, 2) = i;
                    srows(rowsC, 3) = j;
                end
            end
        end
        res = struct('perfect', 0, 'partial80', 0, 'partial50', 0, 'partial0', 0, 'false_positives', 0, 'false_negatives', 0, 'score', 0);
        sortrows(srows(1:rowsC, :), 1);
        
        for i = rowsC:-1:1
            aInd = srows(i, 2);
            bInd = srows(i, 3);
            if (~usedA(aInd) && ~usedB(bInd))
                score = srows(i, 1);
                if (score > 0)
                    usedA(aInd) = true;
                    usedB(bInd) = true;
                    res.score = score + score;
                    if (score == 1)
                        res.perfect = res.perfect + 1;
                    elseif (score >= 0.8)
                        res.partial80 = res.partial80 + 1;
                    elseif (score >= 0.5)
                        res.partial50 = res.partial50 + 1;
                    else
                        res.partial0 = res.partial0 + 1;
                    end
                end
            end
        end
        
        for i = 1:m
            if (~usedA(i))
                res.score = res.score - 1;
                res.false_negatives = res.false_negatives + 1;
            end
        end
        
        for i = 1:n
            if (~usedB(i))
                res.score = res.score - 1;
                res.false_positives = res.false_positives + 1;
            end
        end
    end

    % Returns scores and statistics obtained by genes found by Viterbi.
    % Input:
    %       <allGenes> - annotated genes found by Viterbi
    % Output:
    %       <score>    - structure with score and performance information
    %                    (as returned by match_genes)
    function [score] = get_score(allGenes)
        allCorrectGenes = seq.gene;
        m = length(allCorrectGenes);
        n = length(allGenes);
        genes = [];
        genesC = [];
        correctGenes = [];
        correctCGenes = [];
        
        for i = 1:m
            if (allCorrectGenes(i).Indices(1) < allCorrectGenes(i).Indices(2))
                count = length(correctGenes);
                correctGenes(count + 1).Indices(1) = allCorrectGenes(i).Indices(1);
                correctGenes(count + 1).Indices(2) = allCorrectGenes(i).Indices(2);
            else
                count = length(correctCGenes);
                correctCGenes(count + 1).Indices(1) = allCorrectGenes(i).Indices(1);
                correctCGenes(count + 1).Indices(2) = allCorrectGenes(i).Indices(2);
            end
        end
        
        for i = 1:n
            if (allGenes(i).Indices(1) < allGenes(i).Indices(2))
                count = length(genes);
                genes(count + 1).Indices(1) = allGenes(i).Indices(1);
                genes(count + 1).Indices(2) = allGenes(i).Indices(2);
            else
                count = length(genesC);
                genesC(count + 1).Indices(1) = allGenes(i).Indices(1);
                genesC(count + 1).Indices(2) = allGenes(i).Indices(2);
            end
        end
        
        norm = match_genes(correctGenes, genes);
        comp = match_genes(correctCGenes, genesC);
        score = struct('perfect', norm.perfect + comp.perfect, 'partial80', norm.partial80 + comp.partial80, 'partial50', norm.partial50 + comp.partial50, 'partial0', norm.partial0 + comp.partial0, 'false_positives', norm.false_positives + comp.false_positives, 'false_negatives', norm.false_negatives + comp.false_negatives, 'score', norm.score + comp.score);
    end

    % Simply runs Viterbi and starts matching.
    % Output:
    %       <score> - performance assesment as given by get_score function
    function [score] = run_test()
        fprintf('[i] Running Viterbi algorithm...\n');
        res = hmm.viterbi(seq.Sequence);
        genes = get_genes([], res, false);
        res = hmm.viterbi(seqrcomplement(seq.Sequence));
        genes = get_genes(genes, res, true);
        
        score = get_score(genes);
    end

    % Temporary function used to experiment with parameters of simple
    % intergenic model.
    % Output:
    %       <bestVal> - best value found for transition from simple
    %                   intergenic model back to itself
    function [bestVal] = find_best_value()
        step = 0.01;
        l = eps;
        r = 1 - eps;
        while (step >= 10 ^ -5)
            bestVal = -1;
            bestScore = -Inf;
            for val = l:step:r
                hmm.Trans(1, 1) = val;
                hmm.Trans(1, 2) = 1 - val;
                
                tic;
                testRes = run_test();
                score = testRes.score;
                if (score > bestScore)
                    bestScore = score;
                    bestVal = val;
                end
                
                fprintf('Val = %.5f (score = %.5f)\n', val, score);
                toc;
            end
            l = min([bestVal - step, eps]);
            r = max([bestVal + step, 1 - eps]);
            fprintf('Best value = %.6f\n', bestVal);
            step = step / 10;
        end
    end

    % Making Viterbi work =)
    hmm = HMM(hmm);
    hmm.insert_state(1, 'entry_point', zeros(1, 4));
    hmm.add_edge(1, 'intergenic_dumb', 1); % Change me, that's incorrect...
    
    fprintf('[i] Fixing non-traditional nucleotides in sequence...\n');
    seq.Sequence = fix_sequence(seq.Sequence);
    
    %testRes = find_best_value();
    testRes = run_test();
end