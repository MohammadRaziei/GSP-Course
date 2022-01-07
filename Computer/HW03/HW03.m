clc;clear;close all
addpath("../../gspbox"); gsp_start
addpath("GSP_Computer_HW3");
%%
% n = 12;
% G = gsp_random_regular(n,3);
% G.coords = [cos(2*pi/(n)*(0:n-1)).' sin(2*pi/(n)*(0:n-1)).'];
% gsp_plot_graph(G);
% %%
% [v1,G] = MyVertexSelection(G);
% G1 = MySKReduction(G, v1);
% x = rand(G.N, 1);
% x1 = MyDS(x,v1);
% y = MyInterpolate(G, v1, x1, 0);
%%
table = readtable('Data_City.csv');
isadmin = or(strcmp(table.capital,'primary'), strcmp(table.capital,'admin'));
lat = table.lat(isadmin);
lng = table.lng(isadmin);
n = length(lat);
dist = zeros(n,n);
for i = 1:n
    for j = i+1:n
        dist(i,j) = getDistance(lat(i),lng(i), lat(j),lng(j));
        dist(j,i) = dist(i,j);
    end
end
%%
W = MyDistToGraph(dist,2e5,2e-2);
% full(W)
% length(nonzeros(W))
coords = [lat, lng];
G = gsp_graph(W, coords);
G = gsp_compute_fourier_basis(G);
figure('Color', 'w');
gsp_plot_graph(G)
%%
x = table.Temp(isadmin);
[Gs, ms, ys, xs] = MyPyramidAnalysis(G, x, 3, 0);
xn = xs{end};
xhat = MyPyramidSynthesis(Gs, ms, ys, xn);
figure('Color', 'w');sgtitle("flag=0")
plot_MyPyramidAnalysis(Gs, ys, xs);
figure('Color', 'w');sgtitle("flag=0")
subplot(1,3,1); gsp_plot_signal(G,x); title("orginal signal")
subplot(1,3,2); gsp_plot_signal(G,xhat); title("reconstructed signal")
subplot(1,3,3); gsp_plot_signal(G,x-xhat); title("comparison")

[Gs, ms, ys, xs] = MyPyramidAnalysis(G, x, 3, 1);
xn = xs{end};
xhat = MyPyramidSynthesis(Gs, ms, ys, xn);
figure('Color', 'w');sgtitle("flag=1")
plot_MyPyramidAnalysis(Gs, ys, xs);
figure('Color', 'w');sgtitle("flag=1")
subplot(1,3,1); gsp_plot_signal(G,x); title("orginal signal")
subplot(1,3,2); gsp_plot_signal(G,xhat); title("reconstructed signal")
subplot(1,3,3); gsp_plot_signal(G,x-xhat); title("comparison")

%%
[Gs, ms, ys, xs] = MyPyramidAnalysis(G, x, 3, 0);
xn = awgn(xs{end}, 5, 'measured');
xhat = MyPyramidSynthesis(Gs, ms, ys, xn);
%%
figure('Color', 'w');
subplot(1,3,1); gsp_plot_signal(G,x); title("original signal")
subplot(1,3,2); gsp_plot_signal(G,xhat); title("reconstructed signal")
subplot(1,3,3); gsp_plot_signal(G,abs(x-xhat)); 
title(['E\{||x_{orig} - x_{rec}||^2\} = ' num2str(mean((x-xhat).^2))])
%%
c = 0.01;
xt = @(t) G.U*diag(exp(-c*round(G.e,10)*t))*G.U'*x;
figure('Color', 'w');
for i = 1:5
    subplot(2,3,i); gsp_plot_signal(G,xt(i-1)); title(sprintf("%ith day", i-1))
end
subplot(2,3,6); gsp_plot_signal(G,xt(1e12)); title(sprintf("infinity"))

%%
figure('Color', 'w');
for i = 1:5
    subplot(2,3,i); bar(xt(i-1)); title(sprintf("%ith day", i-1)); ylim([0 20])
end
subplot(2,3,6); bar(xt(1e12)); title(sprintf("infinity")); ylim([0 20])
%%
function plot_MyPyramidAnalysis(Gs, ys, xs)
G0 = Gs{1};
for i = 1:4
    G1 = Gs{i};
    x1 = xs{i};
    subplot(2,4,i);
    gsp_plot_signal(G1, x1); %caxis([0 20]);
    title({sprintf('level %i',i-1) sprintf("x^{(%i)}",i-1)});
    if i ~= 1
        y1 = ys{i-1};
        subplot(2,4,i+4);
        gsp_plot_signal(G0, abs(y1));
%         if(max(abs(y1)) > 20), caxis([0,20]);end
        title({sprintf("y^{(%i)}",i-2), sprintf("||y^{(%i)}||_2^2 = %f",i-2, sum(y1.^2))});
    end
    G0 = G1;
end
end





