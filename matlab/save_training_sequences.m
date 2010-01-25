function [] = save_training_sequences(genome, file)
    seq = get_training_sequences(genome);
    
    f = fopen(file,'w');
    for i=1:length(seq)
        fwrite(f, cell2mat(seq{i}));
        fprintf(f,'\n');
    end
    fclose(f);
end