function Q = Q_per_H(Q,w,j,l)
% Q: matrice
% w: vettore di Householder
% j: colonna corrente
% l: lunghezza vettore w
% 
% L'idea è come per la H_per_A, ma qua c'è un fattore in più da
% considerare, ovvero che un pezzo di Q non è zero e quindi sono due i
% prodotti da fare:
% 
%         | Q_11 Q_12 |   | I    0  |   | Q_11  Q_12*H_t |
% Q * H = |           | * |         | = |                |
%         | Q_21 Q_22 |   | 0   H_t |   | Q_21  Q_22*H_t |
% 


% prime l righe, colonne dalla j-esima alla fine
% Q_12 = Q(1:end-l, j:end); 
% ultime l righe, colonne dalla j-esima alla fine
% Q_22 = Q(end-l+1:end, j:end);

% Q_12*H_t
% Q(1:end-l, j:end) = Q_12-2*(Q_12*w)*w';
% Q(1:l, j:end) = Q_12*H;

% Q_22*H_t
% Q(end-l+1:end, j:end)  = Q_22-2*(Q_22*w)*w';
% Q(end-l+1:end, j:end)  = Q_22*H;

% operazioni riducibili semplicemente a:
Q(:, j:end) = Q(:, j:end) - 2*(Q(:, j:end)*w)*w';

