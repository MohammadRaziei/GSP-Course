function [Gs, ms, ys, xs] = MyPyramidAnalysis(G, x, N, flag)
if nargin < 4, flag = 0; end
Gs = cell(N+1,1); xs = cell(N+1,1); ms = cell(N,1); ys = cell(N,1); 
Gs{1} = G; xs{1} = x;
G1 = G; x1 = x;  
for i = 1:N
    [G1, m, y, x1] = MyAnalysis(G1,x1,flag);
    Gs{i+1} = G1; xs{i+1} = x1;
    ms{i} = m; ys{i} = y;
end
end