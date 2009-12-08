% list ecoli genes that code for proteins (CDS)
function ListProteins(gbf_f)
    for i=1:length(gbf_f.CDS)
        p = gbf_f.CDS(i);
        disp(sprintf('Protein: %s, start: %d, len:%d, protein: %s, function: %s', ...
            p.gene, p.Indices(1), p.Indices(2)-p.Indices(1), p.protein_id, p.function));
    end
end