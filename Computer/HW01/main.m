clc;clear;close all
addpath("../../gspbox"); gsp_start
%% q02
W1 = [0,   0.7, 1.1, 2.3
      0.7, 0,   0,   0
      1.1, 0,   0,   0
      2.3, 0,   0,   0];
coords1 = [0.5,0.41; 0,0; 1,0; 0.5,1];
G1 = gsp_graph(W1, coords1);

W2 = [0,   1.6, 2.4
      1.6, 0,   0.8 
      2.4, 0.8, 0  ];
coords2 = [0.5,0.41; 0,0; 1,0];
G2 = gsp_graph(W2, coords2);

figure("Color", "w")
subplot(121); gsp_plot_graph(G1); title("G1")
subplot(122); gsp_plot_graph(G2); title("G2")


%% q03
param = struct;
param.rule = "cartesian";
G_s = gsp_graph_product(G1,G2, param);
G_s.coords = rand(12,2);
param.rule = "kronecker";
G_t = gsp_graph_product(G1,G2, param);
G_t.coords = rand(12,2);

disp(all(all(full(G_t.A) == kron(full(G1.A), full(G2.A)))))
disp(all(all(G_t.W == kron(G1.W, G2.W))));

disp(all(all(full(G_s.A) == kron(full(G1.A), full(G2.A)))))
disp(all(all(G_s.W == kron(G1.W, G2.W))));

figure("Color", "w")
subplot(121); gsp_plot_graph(G_s); title("cartesian: G_{s}")
subplot(122); gsp_plot_graph(G_t); title("kronecker: G_{t}")

%% q04
myG = G_s;
signal = rand(myG.N, 1);
figure("Color", "w")
gsp_plot_signal(myG,signal); title("myG")
%% q05
myG = gsp_compute_fourier_basis(myG);
% gsp_plot_signal_spectral(G,fhat);
myG.e
myG.U

stem(myG.e); set(gcf, "Color","w")
xlim([0.5 numel(myG.e)+0.5])
xticks(1:numel(myG.e))
%% q06
figure("Color", "w")
i_list = [1 2 11 12];
for c = 1:4
    subplot(2,2,c)
    i = i_list(c);
    signal_i = myG.U(:,i);
    gsp_plot_signal(myG,signal_i); title(sprintf("myG - U_{%i} (e=%f)",c, myG.e(i)));
end

%% part II
%% q07
GL = gsp_logo();
gsp_plot_graph(GL);  
signal = zeros(1,GL.N);
signal(GL.info.idx_g) = -1;
signal(GL.info.idx_s) = 1;
signal(GL.info.idx_p) = -0.5;

figure("Color", "W")
gsp_plot_signal(GL,signal);
%% q08
GL = gsp_compute_fourier_basis(GL);
GL.e;
GL.U;
stem(GL.e)
%% q09
f = GL.U(:,2:3);
figure("Color", "w")
subplot(121); gsp_plot_signal(GL,f(:,1));
subplot(122); gsp_plot_signal(GL,f(:,2));

%% q10
figure("Color", "w")
plot(f(:,1), f(:,2), 'bo')
xlabel("f(:,1)")
ylabel("f(:,2)")


%% q11
rng(1)
idx = kmeans(f,3);
%% q12
colors = [-1, 1, -.5];
sig2 = colors(idx);
figure("Color", "W")
subplot(121); gsp_plot_signal(GL,signal); 
subplot(122); gsp_plot_signal(GL,sig2); title("clustered");
%% q13
rng(1)
f = GL.U(:,2:4);
idx = kmeans(f,3);
colors = [-1, -.5, 1];
sig3 = colors(idx);
figure("Color", "W")
subplot(121); gsp_plot_signal(GL,signal); 
subplot(122); gsp_plot_signal(GL,sig3); title("clustered");