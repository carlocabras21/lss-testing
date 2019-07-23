function [c, s, G] = mygivens(x, i, j)
% i: indice primo numero
% j: indice secondo numero, quello da eliminare

% calcolo di c ed s col metodo "sicuro"
if abs(x(j)) > abs(x(i))
    t = x(i) / x(j);    
    z = sqrt(1 + t*t);
    s = 1/z;
    c = t*s;
else
    t = x(j) / x(i);
    z = sqrt(1 + t*t);
    c = 1/z;
    s = t*c;

end

% costruzione matrice G
if nargout > 2
    G = eye(length(x));

    G(i,i) = c;
    G(i,j) = s;
    
    G(j,i) = -s;
    G(j,j) = c;

end
