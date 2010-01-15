% Downloads genome from GenBank and stores it on disk in a text format
% (one that is a little easier to parse than GBK)
function save_genome(name, filename)
    
    g = get_bank_file(name);
    f = featuresparse(g.Features);
    
    file = fopen(filename, 'w');
    
    write_named_value(g
    
    for i=1:length(f) 
    end
    
    
    fclose(file);

    function write_named_value(v, name)
        fprintf(file,'%s = ', v);
        write_value(v);
        fprintf(file, '\n');
    end
    
    function write_value(v)
        
        if isstruct(v)
            write_struct(v);
        elseif ischar(v)
            fprintf(file, '"%s"', v);
        elseif isnumeric(v) | iscell(v)
            write_array(v);
        elseif isscalar(v)
            fprintf(file,'%f',v); 
        else
            error(['unknown value type: ' v]);
        end
    end

    function write_array(v)
        fprintf(file,'{');
        for i=1:length(v)
            write_value(v(i));
        end
        fprintf(file, '}\n');
    end
   
    function write_struct(s)
        fprintf(file,'{\n');
        names = fieldnames(s);
        for i=1:length(names)
            v = s.(names{i});
            fprintf(file, '%s = ', names{i});
            fprintf('writing %s\n', names{i});
            write_value(v);
            fprintf(file, '\n');
        end
        fprintf(file,'}\n');
    end
    
end