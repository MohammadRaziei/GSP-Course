function x = MyPyramidSynthesis(Gs, ms, ys, xn)
N = length(ys);
x = xn;
for i = N:-1:1
    G = Gs{i}; y = ys{i}; m = ms{i};
    x = MySynthesis(G, m, y, x);
end
end