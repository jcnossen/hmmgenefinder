% download genbank file if not found on disk
function gbk = GetGenbankFile(name)
    file = [name '.gbk'];

    if exist(file, 'file') == 0
        disp(sprintf('Downloading %s from genbank...', name));
        gbk = getgenbank(name, 'TOFILE', file);
    else
        disp('Already downloaded, reading from file...');
        gbk = genbankread(file);
    end
end
