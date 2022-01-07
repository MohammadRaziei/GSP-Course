function L1 = MyKronReductLap(L, v1)
N = length(L);
v1c = setdiff((1:N)', v1);
L1 = L(v1,v1) - L(v1,v1c)*pinv(full(L(v1c, v1c)))*L(v1c,v1);
end