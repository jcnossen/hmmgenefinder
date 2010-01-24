% Downloads genbank file if not found on disk
% Input:
%        name - name of the file
% Output:
%        gbk  - reference to gene bank file.
%
% ------------------------------------------------------------------------
% DBDM - 4, Jelmer Cnossen | Delft University 2009/2010
% ------------------------------------------------------------------------
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
