function R = mychol(A)

% ****** controllo correttezza matrice
% lancia errore se non è simmetrica o se non è definita positiva, ovvero
% se non ha tutti gli autovalori maggiori di zero
if ~issymmetric(A) || ~(all(eig(A)) > 0)
    error("La matrice dev'essere simmetrica definita positiva");
end


% ****** calcolo R
[~, n] = size(A);  % numero elementi
R = zeros(n,n); % pre-allocazione della memoria ed inizializzazione a zero

R(1,1) = sqrt(A(1,1)); % calcolo r_11 fuori dal ciclo

for j = 2:n % dalla seconda colonna in poi
    for i = 1:j-1 % indice i che si ferma prima della diagonale
        % calcolo r_ij, ovvero elementi sopra-diagonale
        R(i,j) = ( A(i,j) - R(1:i-1,i)'*R(1:i-1,j)) / R(i,i);
    end
    
    % calcolo r_jj, ovvero elementi sulla diagonale
    R(j,j) = sqrt( A(j,j) - R(1:j-1,j)'*R(1:j-1,j));
end
