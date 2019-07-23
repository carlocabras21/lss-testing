% TESINA - Computational Mathematics, Carlo Cabras - 60/73/65113


% da testare:
%  - al variare di m ed n
%     - per ogni algoritmo mio e confronto con MATLAB
%        - tempi di calcolo
%        - condizionamento
%        - norma della differenza con la soluzione vera
% 
% 

m = 10;
n = 10;
A = rand(m,n);

% matrice triangolare superiore con due diagonali in più
for j = 1:n
    for i = j+3:m
        A(i,j) = 0;
    end
end

via = 1

tic;
[Q,R] = myqr(A, "householder-light");
tempo1 = toc;

tic;
[Q,R] = myqr(A, "givens-light");
tempo2 = toc;

tic;
[Q,R] = qr(A);
tempo3 = toc;


[tempo1 tempo2 tempo3]