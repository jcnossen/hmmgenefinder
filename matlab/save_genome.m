% Downloads genome from GenBank and stores it on disk in a text format
% (one that is a little easier to parse than GBK)
function save_genome(data, filename)
    if ischar(data)
        g = get_bank_file(name);
        f = featuresparse(g.Features);
    else
        g = data{1};
        f = data{2};
    end
    
    indents = 0;
    
    % write sequence data to a seperate flat text file
    seqfilename = [filename '.seq'];
    seqfile = fopen(seqfilename, 'w');
    fprintf(seqfile, g.Sequence);
    fclose(seqfile);

    file = fopen(filename, 'w');
    write_gbk(g);
    write_features(f);
    fclose(file);
    

    function write_field(sv, name)
        fprintf(file,'%s = ', name);
        write_value(sv.(name));
        fprintf(file, '\n');
    end

    function write_features(f)
        fprintf(file, 'features = ');
        write_struct(f, {'tRNA', 'CDS', 'rRNA', 'misc_RNA', 'misc_feature'});
    end

    function write_gbk(g)
        write_field(g, 'LocusName');
        write_field(g, 'LocusSequenceLength');
        write_field(g, 'LocusTopology');
        write_field(g, 'LocusMoleculeType');
        write_field(g, 'LocusGenBankDivision');
        write_field(g, 'LocusModificationDate');
        write_field(g, 'Definition');
        write_field(g, 'Source');
        tabs(); fprintf(file, 'SequenceFile = %s\n', seqfilename);
    end
    
    function write_value(v)
        [r k]=size(v);
        if k>1
            if iscell(v)
                error('no cell array support implemented');
            elseif r>1
                error('no matrix writing implemented');
            elseif ischar(v)
                fprintf(file, '"%s"', v);
            else
                write_array(v)
            end
        else
            if isstruct(v)
                write_struct(v);
            elseif ischar(v) % single char
                fprintf(file, '"%s"', v);
            elseif isscalar(v)
                fprintf(file,'%f',v); 
            else
                error(['unknown value type: ' v]);
            end
        end
    end

    function write_array(v)
        [rows cols]=size(v);
        if (rows>1)
            error('only row vectors (#rows=1)');
        end
        fprintf(file,'{\n');
        indents=indents+1;
        for i=1:cols
            tabs();
            write_value(v(i));
            fprintf(file,'\n');
        end
        indents=indents-1;
        tabs(); fprintf(file, '}\n');
    end
   
    function write_struct(s, fields)
        fprintf(file,'{\n');
        if nargin == 1
            fields = fieldnames(s);
        end
        indents=indents+1;
        for i=1:length(fields)
            v = s.(fields{i});
            if length(v)==0
                continue
            end
            tabs();
            fprintf(file, '%s = ', fields{i});
            %fprintf('writing %s\n', fields{i});
            write_value(v);
            fprintf(file, '\n');
        end
        indents=indents-1;
        tabs(); fprintf(file,'}\n');
    end
    
    function tabs()
        for i=1:indents
            fprintf(file,'\t');
        end
    end
    
end