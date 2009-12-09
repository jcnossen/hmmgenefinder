% download genbank file if not found on disk
function gbk = get_bank_file(name)
    file = [name '.gbk'];
    
    if (exist(file, 'file') == 0)
        fprintf('[i] Downloading %s from genbank.\n', name);
        gbk = getgenbank(name, 'TOFILE', file);
    else
        fprintf('[i] Sequence %s already exists, reading from file.\n', name);
        gbk = genbankread(file);
    end
end
