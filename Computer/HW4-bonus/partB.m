clc; clear all; close all
addpath("../../gspbox"); gsp_start
addpath("GSP_Computer_HW_Bonus", "GSP_Computer_HW_Bonus/draw")
addpath("MyPyramidTransform")
%%
queries_csv = readtable("queries.csv");
queries_csv.tokens = arrayfun(@(qt) jsondecode(strrep(qt,"'", '"')), queries_csv.tokens,...
    'UniformOutput', false);
num_queries = size(queries_csv,1);
load("tokens_of_queries.mat")
load('queries.mat')
%%
close all
Gqueries = cell(1, num_queries);
for qi = 1:num_queries
    tokens = tokens_of_queries{qi};
    Gqueries{qi} = query2graph(tokens, floor(0.35*size(tokens,1)));
    draw_animal_graph(Gqueries{qi}.L, queries_csv.tokens{qi});
    set(gcf, 'Name', sprintf('query: %g',qi), 'MenuBar', 'none', ...
        'numbertitle', 'off', 'Units', 'normalized', 'Position',[0.4 0.05 0.55 0.9]);
    save_figure(gcf, sprintf('Report/figs/partB-query2graph-%g.png', qi))
end
close all
%%
X = cell(1,num_queries);
for qi = 1:num_queries
    X{qi} = cosine_similarity(tokens_of_queries{qi}, queries{qi});
    VV = GCModulMax1(Gqueries{qi}.W);
    VVmax = max(VV);
    Vmask = VV == unique(VV)';
    Vmax = max(sum(Vmask,1));
    [~,I] = max(X{qi}.*Vmask,[],1);
    tokens = queries_csv.tokens{qi};
    clusters = cell(Vmax+1, VVmax);
    for v = 1:VVmax
        mask = Vmask(:,v);
        clusters(:,v) = [tokens(I(v)); tokens(mask); repmat({' '}, Vmax - sum(mask),1)];
    end
    writetable(cell2table(clusters(2:end,:),'VariableNames',clusters(1,:)), ...
        sprintf('Report/figs/partB-clusters-%g.csv',qi))
end
%%
qi = 1;
[Gs, ms, ys, xs] = MyPyramidAnalysis(Gqueries{qi}, X{qi}, 1);







