function y = MyInterpolate(G, v1, x, flag, epsilon)
if nargin < 4, flag = 0; end
if nargin < 5, epsilon = 0.005; end
if flag == 0
    Phi = MyInterpBasis(G, v1, epsilon);
    alpha = Phi(v1,:)\x;
    y = Phi*alpha;
else
    k = length(v1);
    y = G.U(:,1:k)*(G.U(v1,1:k)\x);
end
end