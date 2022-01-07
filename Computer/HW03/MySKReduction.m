function G1 = MySKReduction(G, v1, Q)
G1 = MyKronReduction(G, v1);
if nargin < 3, Q = round(4*G1.N*log2(G1.N)); end
G1 = MySpectralSparsification(G1, Q);
end