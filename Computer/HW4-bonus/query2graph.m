function Gquery = query2graph(tokens, k, normalized)
if nargin < 2, k = 5; end
if nargin < 3, normalized = true; end
CS = cosine_similarity(tokens);
Z = zero_diag(1 - CS).^2;% cosine distance
% Z = (1 ./ (abs(CS) + 1e-8) )-1;
theta = gsp_compute_graph_learning_theta(Z, k);
W = gsp_learn_graph_log_degrees(theta*Z, 1, 1);
W(W<1e-3) = 0; % clean up zeros
if normalized && max(W(:)) > 1e-3, W = W / max(W(:)); end
G = gsp_ring(size(W,1));
Gquery = gsp_graph(W, G.coords);
end


