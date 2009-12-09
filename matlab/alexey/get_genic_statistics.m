function [sum_arr, sum_nuc] = get_genic_statistics(seq, traditional_only)
    if (nargin < 2)
        traditional_only = false;
    end

    n = length(seq.gene);
    sum_arr = zeros([4 4 4]);
    sum_nuc = zeros(1, 16);
    for i = 1:n
        s = min(seq.gene(i).Indices);
        f = max(seq.gene(i).Indices);
        str = seq.Sequence(s:f);
        int_str = nt2int(str);
        if (~traditional_only || sum(int_str > 4) == 0)
            [tmp, tmp_arr] = codoncount(str, 'Reverse', seq.gene(i).Indices(1) > seq.gene(i).Indices(2));
            sum_arr = sum_arr + tmp_arr;
            sum_nuc(int_str) = sum_nuc(int_str) + 1;
        end
    end
end