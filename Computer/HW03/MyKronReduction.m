function G1 = MyKronReduction(G, v1)
L1 = MyKronReductLap(G.L, v1);
W1 = diag(diag(L1)) - L1;
coords = G.coords(v1,:);
G1 = gsp_graph(W1, coords);
end