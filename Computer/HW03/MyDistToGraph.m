function W = MyDistToGraph(dist, tau, thresh)
    W = exp(-(dist.^2)/tau);
    W = W - diag(diag(W));
    W = sparse(W.*double(W>thresh));
end