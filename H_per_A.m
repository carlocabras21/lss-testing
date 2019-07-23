function A = H_per_A(A,w,j,l)
% A: matrice
% w: vettore di Householder
% j: colonna corrente
% l: lunghezza vettore w
%
%         | I    0  |   | A_11 A_12 |   | A_11    A_12   |
% H * A = |         | * |           | = |                |
%         | 0   H_t |   |  0   A_22 |   |   0   A_22*H_t |
%
% questa funzione implementa H*A, quando la H è stata creata "bordando"
% la H_tilde. La "bordatura" è implementata andando a modificare soltanto
% soltanto la parte in basso a destra della matrice. Non ci sono prodotti
% tra matrici (compl. n^3) ma soltanto prodotti tra vettori, quindi 
% massimo n^2 cioè quando si ha un prodotto colonna*riga o matrice*vettore.

% A_22 è la matrice con le ultime l righe e n-j colonne
A_22 = A(end-l+1:end, j:end);

% modifico solo quelle parti
A(end-l+1:end, j:end)  = A_22-2*w*(w'*A_22);
