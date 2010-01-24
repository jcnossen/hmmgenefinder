function [] = save_training_sequences(genome, file)
    seq = HMM.get_training_sequences(genome);
    
    f = fopen(file,'w');
    for i=1:length(seq)
        fwrite(f, seq{i});
        if i==length(seq)
            fwrite(f, ',\n');
        end
    end
    fclose(f);
end