function G1 = MySpectralSparsification(G,Q)
[I,J,We] = find(G.W);
n = length(We);
P = zeros(n,1);
for k = 1:n
    P(k) = MyResistanceDistance(G,I(k),J(k))*We(k);
end
P = P / sum(P);
PS = cumsum(P);
% I1 = zeros(Q,1); J1 = zeros(Q,1);
% We1 = zeros(Q,1);
W1 = zeros(G.N, G.N);
for q = 1:Q
    s = 1+sum(PS < rand);
%     I1(q)= I(s);
%     J1(q)= J(s);
%     We1(q)= We(s);
    W1(I(s), J(s)) = W1(I(s), J(s)) + G.W(I(s), J(s))/(Q*P(s));
end
% W1 = sparse(I1,J1,We1,G.N,G.N);
G1 = gsp_graph(W1, G.coords);
end