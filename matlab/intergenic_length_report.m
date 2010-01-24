% Outputs some statistics on lengths of intergenic regions
% Input:
%       <seq> - sequence with intergenic regions annotated
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
function [] = intergenic_length_report(seq)
    seq = HMM_Intergenic.get_intergenic(seq, false, false);
    
    n = length(seq.gene);
    lenA = zeros(1, n);
    for i = 1:n
        len = abs(seq.gene(i).Indices(1) -  seq.gene(i).Indices(2)) + 1;
        lenA(i) = len;
    end
    fprintf('Min length: %i\n', min(lenA));
    fprintf('Max length: %i\n', max(lenA));
    fprintf('\n');
    fprintf('Total intergenic: %i\n', n);
    fprintf('Length 1 - 14: %i\n', sum((lenA <= 14) .* (lenA >= 1)));
    fprintf('Length 10 - Inf: %i\n', sum(lenA >= 10));
end