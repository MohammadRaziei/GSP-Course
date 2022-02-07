clc; clear all; close all
addpath("../../gspbox"); gsp_start
addpath("GSP_Computer_HW_Bonus", "MyPyramidTransform")
%%
queries_csv = readtable("queries.csv");
queries_csv.tokens = arrayfun(@(qt) jsondecode(strrep(qt,"'", '"')), queries_csv.tokens,...
    'UniformOutput', false);
num_queries = size(queries_csv,1);
load("tokens_of_queries.mat")
load('queries.mat')
load('docs.mat'); docs = cell2mat(docs);
%%
Gdocs = query2graph(docs, floor(0.4*size(docs,1)));
Gdocs.Ne % 38761
figure('Color','w')
imagesc(Gdocs.W); colorbar
save_figure(gcf, 'Report/figs/partC-Gdocs-W.png');
%%
% param = struct();
% param.show_edges = 1;
% gsp_plot_graph(Gdocs, param)
%%
X = cosine_similarity(docs, queries{1});
[Gs, ms, ys, xs] = MyPyramidAnalysis(Gdocs, X, 7);
Gselected = Gs{end};
N = Gselected.N; % 17
Xselected = xs{end};
Mselected = ms{end};
Xsampled = X(Mselected);
disp((Xsampled/max(X))') 
figure('Color','w')
gsp_plot_signal(Gdocs, X);
save_figure(gcf, 'Report/figs/partC-Gdocs.png');
figure('Color','w')
gsp_plot_signal(Gselected, Xselected);
save_figure(gcf, 'Report/figs/partC-Gselected.png');

%%
Gsampled = query2graph(docs(Mselected,:), floor(0.4*N), false);
CS = cosine_similarity(docs(Mselected,:));
norm(1-CS, 'fro')
var(CS,0,'all')
draw_animal_graph(Gsampled.L, arrayfun(@num2str, find(Mselected), 'UniformOutput', 0))
set(gcf, 'MenuBar', 'none', 'numbertitle', 'off', 'Units', 'normalized', ...
    'Position',[0.4 0.05 0.55 0.9]);
save_figure(gcf, 'Report/figs/partC-Gsampled.png');
