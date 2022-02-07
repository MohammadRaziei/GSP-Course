function cs = cosine_similarity(x,y)%row data
if nargin < 2, y = x; end
n = size(x, 1);
m = size(y, 1);
normX = sqrt(sum(abs(x).^2, 2));
normY = sqrt(sum(abs(y).^2, 2));
cs = zeros(n,m);
for i = 1:n
   for j = 1:m
       cs(i,j) = x(i,:)*y(j,:)' / (normX(i) * normY(j));
   end
end
end