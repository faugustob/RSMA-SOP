function code = classify4(a, b, c, d)
%CLASSIFY4 Classifies 4 numbers by equality pattern
% Returns:
%   4    -> all four equal
%   31   -> three equal, one different
%   22   -> two distinct pairs
%   211  -> one pair, two distinct singles
%   1111 -> all four different

    x = [a b c d];
    counts = sort(histcounts(categorical(x)), 'descend');

    if isequal(counts, [4])
        code = 4;
    elseif isequal(counts, [3 1])
        code = 31;
    elseif isequal(counts, [2 2])
        code = 22;
    elseif isequal(counts, [2 1 1])
        code = 211;
    else
        code = 1111;
    end
end
