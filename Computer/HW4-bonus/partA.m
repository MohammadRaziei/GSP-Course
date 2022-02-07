clc; clear all; close all
addpath("../../gspbox"); gsp_start
addpath("GSP_Computer_HW_Bonus", "GSP_Computer_HW_Bonus/draw")
%%
load('Data.mat')
Z = gsp_distanz(data').^2;
%%
a = 5; b = 5;
[W1] = gsp_learn_graph_log_degrees(Z, a, b);
W1(W1<1e-3) = 0; % clean up zeros
G1 = gsp_graph(W1);
draw_animal_graph(G1.L,names);
set(gcf, 'Units', 'normalized','NumberTitle', 'off', ...
    'outerposition',[0.45 0.1 0.53 .9]);
save_figure(gcf, sprintf('Report/figs/test-gsp-learn-graph-log-degrees-Z-%g-%g.png', a, b));
%%
a = 5;
[W2] = gsp_learn_graph_l2_degrees(Z, a);
W2(W2<1e-3) = 0; % clean up zeros
G2 = gsp_graph(W2);
draw_animal_graph(G2.L,names);
set(gcf, 'Units', 'normalized','NumberTitle', 'off', ...
    'outerposition',[0.45 0.1 0.53 .9]);
save_figure(gcf, sprintf('Report/figs/test-gsp-learn-graph-12-degrees-Z-%g.png', a));
%%
close all
a = 5; b = 5;
[W] = gsp_learn_graph_log_degrees(Z, a, b);
if max(max(W)) > 1e-3,  W = W / max(max(W)); end
W(W<1e-3) = 0; % clean up zeros
G = gsp_update_weights(G1,W);
draw_animal_graph(G.L,names);
text(0.3,-0.05,sprintf('gsp-learn-graph-log-degrees(Z, %g, %g)', a, b), ...
    "Units", "normalized")
set(gcf, 'Units', 'normalized','NumberTitle', 'off', ...
    'outerposition',[0.45 0.1 0.53 .9]);
save_figure(gcf, sprintf('Report/figs/gsp-learn-graph-log-degrees-Z-%g-%g.png', a, b));
%%
k = 4;
theta = gsp_compute_graph_learning_theta(Z, k); 
[W] = gsp_learn_graph_log_degrees(theta * Z, 1, 1);
if max(max(W)) > 1e-3,  W = W / max(max(W)); end
W(W<1e-3) = 0; % clean up zeros
G = gsp_update_weights(G1,W);
draw_animal_graph(G.L,names);
text(0.3,-0.05,sprintf('gsp-learn-graph-log-degrees(%0.3g Z,1,1)', theta), ...
    "Units", "normalized")
set(gcf, 'Units', 'normalized','NumberTitle', 'off', ...
    'outerposition',[0.45 0.1 0.53 .9]);
save_figure(gcf, sprintf('Report/figs/gsp-learn-graph-log-degrees-Z-k%g.png', k));

