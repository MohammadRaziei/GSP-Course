function ms = MyFixIndices(ms)
n = size(ms{1},1);
ms{1} = logical(ms{1});
for i = 2:size(ms,1)
    m = zeros(n,1);
    m(ms{i-1}) = ms{i}; 
    ms{i} = logical(m);
end
end