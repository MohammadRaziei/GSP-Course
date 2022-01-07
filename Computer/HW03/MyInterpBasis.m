function Phi = MyInterpBasis(G,v1,epsilon)
if nargin < 3, epsilon = 0.005; end
Phi = G.U*diag(1./(epsilon+G.e))*(G.U(v1,:))';
end