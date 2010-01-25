% Loads E.Coli genome from gene bank and extracts necessary information
% from it (sequence and annotated genes).
% Output:
%       <res> - annotated sequence in format used everywhere in out scripts
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
function [res] = load_ecoli()
    g = get_bank_file('AE005174');
    f = featuresparse(g.Features);
    res.Sequence = g.Sequence;
    n = length(f.gene);
    for i = 1:n
        res.gene(i).Indices = f.gene(i).Indices;
    end
end