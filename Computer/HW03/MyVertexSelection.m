function [vbig, G] = MyVertexSelection(G)
if not(isfield(G, 'U')), G = gsp_compute_fourier_basis(G); end
v1 = find(G.U(:,end) >= 0);
v1c = find(G.U(:,end) < 0);
if numel(v1) >= numel(v1c)
    vbig = v1;
else
    vbig = v1c;
end
end