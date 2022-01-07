function d = MyResistanceDistance(G,i,j)
delta = zeros(G.N,1);
delta([i,j]) = [1 -1];
d = delta'*pinv(G.L)*delta;
end