function x = MySynthesis(G, m, y, x1)
v1 = find(m);
xhat = MyInterpolate(G, v1, x1, 0);
x = y + xhat;
end