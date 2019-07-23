function [w,k,H] = housemat(x)

n = length(x);

sigma = norm(x);
if x(1) >=0
    k = -sigma;
else
    k = sigma;
end
lambda = sqrt(2 * sigma * (sigma + abs(x(1))));
w = x;
w(1) = w(1) - k;
w = w / lambda;

if nargout > 2
    H = eye(n) - 2*w*w';
end



