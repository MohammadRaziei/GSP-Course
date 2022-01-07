function [y, G] = MyHfilter(G,x)
if not(isfield(G, 'U')), G = gsp_compute_fourier_basis(G); end
y = G.U*diag(1./(1+2*G.e))*G.U'*x;
end