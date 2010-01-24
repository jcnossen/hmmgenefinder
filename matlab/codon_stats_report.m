% Calculates statistics for codon usage in given contig.
% Input:
%       seq - sequence with annotated genes
% Output:
%       res - statistics in form of a structure
% ------------------------------------------------------------------------
% DBDM - 4, Alexey Gritsenko | Leiden University 2009/2010
% ------------------------------------------------------------------------
function [str] = codon_stats_report(seq)
    
    function [res] = output2csv()
        directTotal = zeros(1, 3);
        complementTotal = zeros(1, 3);
        
        endl = sprintf('\n');
        res = [';RF;Train', endl];
        res = [res, 'Codon;;1;2;3;4;5;6', endl];
        fields = fieldnames(direct(1));
        nFields = length(fields);
        
        for codonI = 1:nFields
            codonName = fields{codonI};
            for frameCSV = 1:3
                directTotal(frameCSV) = directTotal(frameCSV) + direct(frameCSV).(codonName);
                complementTotal(frameCSV) = complementTotal(frameCSV) + complement(frameCSV).(codonName);
            end
        end
        
        for codonI = 1:nFields
            codonName = fields{codonI};
            res = [res, codonName, ';;'];
            for frameCSV = 1:3
                res = [res, int2str(direct(frameCSV).(codonName)), ',', sprintf('%.4f', direct(frameCSV).(codonName) / directTotal(frameCSV)), ';'];
            end
            for frameCSV = 1:3
                res = [res, int2str(complement(frameCSV).(codonName)), ',', sprintf('%.4f', complement(frameCSV).(codonName) / complementTotal(frameCSV)), ';'];
            end
            res = [res, endl];
        end
        res = [res, 'Total;;'];
        for frameCSV = 1:3
            res = [res, int2str(directTotal(frameCSV)), ';'];
        end
        for frameCSV = 1:3
            res = [res, int2str(complementTotal(frameCSV)), ';'];
        end
    end
    
    str = seq.Sequence;
    for frame = 1:3
        direct(frame) = codoncount(str, 'Frame', frame);
    end
    
    str = seqrcomplement(seq.Sequence);
    for frame = 1:3
        complement(frame) = codoncount(str, 'Frame', frame);
    end
    str = output2csv;
end