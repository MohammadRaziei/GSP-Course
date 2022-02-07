function [G1, m, y, x1] = MyAnalysis(G,x,flag)
if nargin < 3, flag = 0; end
[v1,G] = MyVertexSelection(G);
m = zeros(G.N, 1); m(v1) = 1;
[x1, G] = MyHfilter(G, x);
x1 = MyDS(x1,v1);
G1 = MyKronReduction(G, v1);
G1 = gsp_compute_fourier_basis(G1);
xhat = MyInterpolate(G, v1, x1, flag);
y = x - xhat;
end