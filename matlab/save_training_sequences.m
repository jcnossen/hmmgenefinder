% Saves training sequences for GHMM
% Input:
%       <seq>  - cell array with training data
%       <file> - file name to save data to
% ------------------------------------------------------------------------
% DBDM - 4, Jelmer Cnossen | Delft University 2009/2010
% ------------------------------------------------------------------------
function [] = save_training_sequences(seq, file)    
    f = fopen(file,'w');
    for i = 1:length(seq)
        fwrite(f, cell2mat(seq{i}));
        fprintf(f,'\n');
    end
    fclose(f);
end