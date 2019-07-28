% TESINA - Computational Mathematics, Carlo Cabras - 60/73/65113

% IMPORTANTE: lanciare i test uno alla volta 

% MAPPA DEI COLORI DELLE LINEE DEI GRAFICI
% blu: Householder (tratteggiato: light)
% rosso: Givens (tratteggiato: light)
% magenta: Cholesky  (tratteggiato: di Matlab)
% verde: qr Matlab (tratteggiato: backslash, puntini: pseudoinversa)

% ********************************************************************** %
% ********************************************************************** %
% ***** PARTE 1: sistemi quadrati
% ********************************************************************** %
% ********************************************************************** %

% ********************************************************************** %
% TEST 1.1: confronto tra QR mio (fattorizzazione di householder e di
% givens, con le loro versioni "light"), QR di Matlab e backslash.
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 50:10:100;

t_sol  = zeros(length(dimensions),6);
e_sol  = zeros(length(dimensions),6);
e_fact = zeros(length(dimensions),5);
e_orth = zeros(length(dimensions),5);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(n,n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    norm_sol = norm(sol);
    b = A*sol;
    
    tic;
    [Q,R] = myqr(A, "householder");
    y = Q'*b;
    x = R\y;
    t_sol(i,1)  = toc;
    e_sol(i,1)  = norm(x-sol)/norm_sol;
    e_orth(i,1) = norm(Q*Q'-eye(n));
    e_fact(i,1) = norm(A-Q*R);
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,2) = toc;
    e_sol(i,2) = norm(x-sol)/norm_sol;
    e_orth(i,2) = norm(Q*Q'-eye(n));
    e_fact(i,2) = norm(A-Q*R);
    
    tic;
    [Q,R] = myqr(A, "givens");
    y = Q'*b;
    x = R\y;
    t_sol(i,3) = toc;
    e_sol(i,3) = norm(x-sol)/norm_sol;
    e_orth(i,3) = norm(Q*Q'-eye(n));
    e_fact(i,3) = norm(A-Q*R);
    
    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,4) = toc;
    e_sol(i,4) = norm(x-sol)/norm_sol;
    e_orth(i,4) = norm(Q*Q'-eye(n));
    e_fact(i,4) = norm(A-Q*R);
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,5) = toc;
    e_sol(i,5) = norm(x-sol)/norm_sol;
    e_orth(i,5) = norm(Q*Q'-eye(n));
    e_fact(i,5) = norm(A-Q*R);
    
    tic;
    x = A\b;
    t_sol(i,6) = toc;
    e_sol(i,6) = norm(x-sol)/norm_sol;
    
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol(:,1), 'b', ...
         dimensions, e_sol(:,2), 'b--', ...
         dimensions, e_sol(:,3), 'r', ...
         dimensions, e_sol(:,4), 'r--', ...
         dimensions, e_sol(:,5), 'g', ...
         dimensions, e_sol(:,6), 'g--');
title("errori di soluzione relativi");
legend('householder','householder-light','givens','givens-light','qr','backslash');

figure(2); 
semilogy(dimensions, t_sol(:,1), 'b', ...
         dimensions, t_sol(:,2), 'b--', ...
         dimensions, t_sol(:,3), 'r', ...
         dimensions, t_sol(:,4), 'r--', ...
         dimensions, t_sol(:,5), 'g', ...
         dimensions, t_sol(:,6), 'g--');
title("tempi di risoluzione");
legend('householder','householder-light','givens','givens-light','qr','backslash');


figure(3); 
semilogy(dimensions, e_orth(:,1), 'b', ...
         dimensions, e_orth(:,2), 'b--', ...
         dimensions, e_orth(:,3), 'r', ...
         dimensions, e_orth(:,4), 'r--', ...
         dimensions, e_orth(:,5), 'g');
title("errori di ortogonalizzazione");
legend('householder','householder-light','givens','givens-light','qr');

figure(4); 
semilogy(dimensions, e_fact(:,1), 'b', ...
         dimensions, e_fact(:,2), 'b--', ...
         dimensions, e_fact(:,3), 'r', ...
         dimensions, e_fact(:,4), 'r--', ...
         dimensions, e_fact(:,5), 'g');
title("errori di fattorizzazione");
legend('householder','householder-light','givens','givens-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");


% ********************************************************************** %
% TEST 1.2: Visti i tempi più alti per Householder e Givens, uso solo le
% loro versioni light. Così posso usare n più grandi.
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 50:50:300;

t_sol  = zeros(length(dimensions),4);
e_sol  = zeros(length(dimensions),4);
e_orth = zeros(length(dimensions),3);
e_fact = zeros(length(dimensions),3);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(n,n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    norm_sol = norm(sol);
    b = A*sol;
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,1) = toc;
    e_sol(i,1) = norm(x-sol)/norm_sol;
    e_orth(i,1) = norm(Q*Q'-eye(n));
    e_fact(i,1) = norm(A-Q*R);

    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,2) = toc;
    e_sol(i,2) = norm(x-sol)/norm_sol;
    e_orth(i,2) = norm(Q*Q'-eye(n));
    e_fact(i,2) = norm(A-Q*R);
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,3) = toc;
    e_sol(i,3) = norm(x-sol)/norm_sol;
    e_orth(i,3) = norm(Q*Q'-eye(n));
    e_fact(i,3) = norm(A-Q*R);
    
    tic;
    x = A\b;
    t_sol(i,4) = toc;
    e_sol(i,4) = norm(x-sol)/norm_sol;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol(:,1), 'b--', ...
         dimensions, e_sol(:,2), 'r--', ...
         dimensions, e_sol(:,3), 'g', ...
         dimensions, e_sol(:,4), 'g--');
title("errori di soluzione relativi");
legend('householder-light','givens-light','qr','backslash');

figure(2); 
semilogy(dimensions, t_sol(:,1), 'b--', ...
         dimensions, t_sol(:,2), 'r--', ...
         dimensions, t_sol(:,3), 'g', ...
         dimensions, t_sol(:,4), 'g--');
title("tempi di risoluzione");
legend('householder-light','givens-light','qr','backslash');

figure(3); 
semilogy(dimensions, e_orth(:,1), 'b--', ...
         dimensions, e_orth(:,2), 'r--', ...
         dimensions, e_orth(:,3), 'g');
title("errori di ortogonalizzazione");
legend('householder-light','givens-light','qr');

figure(4); 
semilogy(dimensions, e_fact(:,1), 'b--', ...
         dimensions, e_fact(:,2), 'r--', ...
         dimensions, e_fact(:,3), 'g');
title("errori di fattorizzazione");
legend('householder-light','givens-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");


% ********************************************************************** %
% TEST 1.3: comparazione tra mio Householder-light e qr di MATLAB, con
% n grandi. A matrice quadrata nxn
% ********************************************************************** %

% * * * * DA DECIDERE SE NECESSARIO * * * * 
% anzi, decidere se usare direttamente questo e bogare il precedente. Mi sa
% che farò così.

% clear variables;
% 
% i = 1;
% 
% dimensions = 1000:500:1500;
% 
% t_sol  = zeros(length(dimensions),2);
% e_sol  = zeros(length(dimensions),2);
% e_orth = zeros(length(dimensions),2);
% e_fact = zeros(length(dimensions),2);
% n_cond = zeros(length(dimensions),1);
% 
% for n = dimensions
%     fprintf('n: %d\n',n);
%     
% 	A = rand(n,n);
%     n_cond(i) = cond(A);
%     
%     sol = ones(n,1);
%     norm_sol = norm(sol);
%     b = A*sol;
%     
%     k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
%            % tenere conto di quanti algoritmi sto testando
%     
%     tic;
%     [Q,R] = myqr(A, "householder-light");
%     y = Q'*b;
%     x = R\y;
%     t_sol(i,k) = toc;
%     e_sol(i,k) = norm(x-sol)/norm_sol;
%     e_orth(i,k) = norm(Q*Q'-eye(n));
%     e_fact(i,k) = norm(A-Q*R);
%     k = k+1;
% 
%     tic;
%     [Q,R] = qr(A);
%     y = Q'*b;
%     x = R\y;
%     t_sol(i,k) = toc;
%     e_sol(i,k) = norm(x-sol)/norm_sol;
%     e_orth(i,k) = norm(Q*Q'-eye(n));
%     e_fact(i,k) = norm(A-Q*R);
%     k = k + 1;
%     
%     i = i+1;
% 
% end
% 
% figure(1); 
% semilogy(dimensions, e_sol); title("errori di soluzione relativi");
% legend('householder-light','qr');
% 
% figure(2); 
% semilogy(dimensions, t_sol); title("tempi di risoluzione");
% legend('householder-light','qr');
% 
% figure(3); 
% semilogy(dimensions, e_orth); title("errori di ortogonalizzazione");
% legend('householder-light','qr');
% 
% figure(4); 
% semilogy(dimensions, e_fact); title("errori di fattorizzazione");
% legend('householder-light','qr');
% 
% figure(5); 
% plot(dimensions, n_cond); title("numero di condizionamento");



% ********************************************************************** %
% TEST 1.4: matrice triangolare superiore con diagonali in più
% dimensioni di A fissate a 500x500, cambia il numero delle sottodiagonali.
% Confrontro tra Householder-light, Givens-light e qr di MATLAB.
% ********************************************************************** %

% mi aspetto che givens-light funziona meglio, visto che se un elemento è
% già zero non calcolo nulla.

% funziona male per matrici quasi triangolari, perché mi dà l'errore:
% "Warning: Matrix is close to singular or badly scaled. Results may be 
% inaccurate. RCOND =  6.550039e-22." quindi parto con 4 sottodiagonali.
% Si nota comunque che Givens-light vince su Householder-light, ma non su
% qr di Matlab

clear variables;

i = 1;

dimensions = 4:2:20;

t_sol  = zeros(length(dimensions),4);
e_sol  = zeros(length(dimensions),4);
e_orth = zeros(length(dimensions),3);
e_fact = zeros(length(dimensions),3);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('sottodiagonali: %d\n',n);
    
    dim = 300;
	A = rand(dim,dim);
    for j = 1:dim
        for k = j+n+1:dim
            A(k,j) = 0;
        end
    end
    n_cond(i) = cond(A);
    
    sol = ones(dim,1);
    norm_sol = norm(sol);
    b = A*sol;

    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,1) = toc;
    e_sol(i,1) = norm(x-sol)/norm_sol;
    e_orth(i,1) = norm(Q*Q'-eye(dim));
    e_fact(i,1) = norm(A-Q*R);
    
    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,2) = toc;
    e_sol(i,2) = norm(x-sol)/norm_sol;
    e_orth(i,2) = norm(Q*Q'-eye(dim));
    e_fact(i,2) = norm(A-Q*R);
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,3) = toc;
    e_sol(i,3) = norm(x-sol)/norm_sol;
    e_orth(i,3) = norm(Q*Q'-eye(dim));
    e_fact(i,3) = norm(A-Q*R);
    
    tic;
    x = A\b;
    t_sol(i,4) = toc;
    e_sol(i,4) = norm(x-sol)/norm_sol;
    
    i = i+1;

end
figure(1); 
semilogy(dimensions, e_sol(:,1), 'b--', ...
         dimensions, e_sol(:,2), 'r--', ...
         dimensions, e_sol(:,3), 'g', ...
         dimensions, e_sol(:,4), 'g--');
title("errori di soluzione relativi");
legend('householder-light','givens-light','qr','backslash');

figure(2); 
semilogy(dimensions, t_sol(:,1), 'b--', ...
         dimensions, t_sol(:,2), 'r--', ...
         dimensions, t_sol(:,3), 'g', ...
         dimensions, t_sol(:,4), 'g--');
title("tempi di risoluzione");
legend('householder-light','givens-light','qr','backslash');

figure(3); 
semilogy(dimensions, e_orth(:,1), 'b--', ...
         dimensions, e_orth(:,2), 'r--', ...
         dimensions, e_orth(:,3), 'g');
title("errori di ortogonalizzazione");
legend('householder-light','givens-light','qr');

figure(4); 
semilogy(dimensions, e_fact(:,1), 'b--', ...
         dimensions, e_fact(:,2), 'r--', ...
         dimensions, e_fact(:,3), 'g');
title("errori di fattorizzazione");
legend('householder-light','givens-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");

% ********************************************************************** %
% ********************************************************************** %
% ***** PARTE 2: sistemi sovradimensionati
% ********************************************************************** %
% ********************************************************************** %

% ********************************************************************** %
% TEST 2.1: Equazioni normali. A m x n con m > n, sitema Ax = b impostato 
% come A'Ax = A'b.
% differenza tra Cholesky mio, chol di Matlab, pseudoinversa e \.
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 10:10:100;

t_sol = zeros(length(dimensions),4);
e_sol = zeros(length(dimensions),4);
e_fact = zeros(length(dimensions),2);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(110,n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    norm_sol = norm(sol);
    b = A*sol;

    tic;
    R = mychol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    t_sol(i,1) = toc;
    e_sol(i,1) = norm(x-sol)/norm_sol;
    e_fact(i,1) = norm(A'*A-R'*R);
    
    tic;
    R = chol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    t_sol(i,2) = toc;
    e_sol(i,2) = norm(x-sol)/norm_sol;
    e_fact(i,2) = norm(A'*A-R'*R);
    
    tic;
    x = pinv(A)*b;
    t_sol(i,3) = toc;
    e_sol(i,3) = norm(x-sol)/norm_sol;
    
    tic;
    x = A\b;
    t_sol(i,4) = toc;
    e_sol(i,4) = norm(x-sol)/norm_sol;
    
    i = i+1;
end

% mychol, chol, pinv, backslash

figure(1); 
semilogy(dimensions, e_sol(:,1), 'm', ...
         dimensions, e_sol(:,2), 'm--', ...
         dimensions, e_sol(:,3), 'g:', ...
         dimensions, e_sol(:,4), 'g--');
title("errori di soluzione relativi");
legend('mychol','chol','pinv','backslash');

figure(2); 
semilogy(dimensions, t_sol(:,1), 'm', ...
         dimensions, t_sol(:,2), 'm--', ...
         dimensions, t_sol(:,3), 'g:', ...
         dimensions, t_sol(:,4), 'g--');
title("tempi di risoluzione");
legend('mychol','chol','pinv','backslash');

figure(4); 
semilogy(dimensions, e_fact(:,1), 'm', ...
         dimensions, e_fact(:,2), 'm--');
title("errori di fattorizzazione");
legend('mychol','chol');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");


% ********************************************************************** %
% TEST 2.2: differenza tra myqr e qr di MATLAB. A matrice rettangolare 
% m x n con m > n.
% Uso la mia versione di qr con householder-light, comparandola alla qr di
% Matlab, alla pseudoinversa ed al backslash.
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 50:50:200;

t_sol  = zeros(length(dimensions),4);
e_sol  = zeros(length(dimensions),4);
e_orth = zeros(length(dimensions),2);
e_fact = zeros(length(dimensions),2);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
    m = 250;
	A = rand(m,n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    norm_sol = norm(sol);
    b = A*sol;
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,1) = toc;
    e_sol(i,1) = norm(x-sol)/norm_sol;
    e_orth(i,1) = norm(Q*Q'-eye(m));
    e_fact(i,1) = norm(A-Q*R);

    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,2) = toc;
    e_sol(i,2) = norm(x-sol)/norm_sol;
    e_orth(i,2) = norm(Q*Q'-eye(m));
    e_fact(i,2) = norm(A-Q*R);
        
    tic;
    x = pinv(A)*b;
    t_sol(i,3) = toc;
    e_sol(i,3) = norm(x-sol)/norm_sol;
        
    tic;
    x = A\b;
    t_sol(i,4) = toc;
    e_sol(i,4) = norm(x-sol)/norm_sol;
        
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol(:,1), 'b--', ...
         dimensions, e_sol(:,2), 'g', ...
         dimensions, e_sol(:,3), 'g:', ...
         dimensions, e_sol(:,4), 'g--');
title("errori di soluzione relativi");
legend('householder-light','qr','pinv','backslash');

figure(2); 
semilogy(dimensions, t_sol(:,1), 'b--', ...
         dimensions, t_sol(:,2), 'g', ...
         dimensions, t_sol(:,3), 'g:', ...
         dimensions, t_sol(:,4), 'g--');
title("tempi di risoluzione");
legend('householder-light','qr','pinv','backslash');

figure(3); 
semilogy(dimensions, e_orth(:,1), 'b--', ...
         dimensions, e_orth(:,2), 'g');
title("errori di ortogonalizzazione");
legend('householder-light','qr');

figure(4); 
semilogy(dimensions, e_fact(:,1), 'b--', ...
         dimensions, e_fact(:,2), 'g');
title("errori di fattorizzazione");
legend('householder-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");

% ********************************************************************** %
% TEST 2.3: matrice triangolare superiore con diagonali in più. Le
% dimensioni di A  sono fissate a 800x300, cambia il numero delle 
% sottodiagonali. Confrontro tra Householder-light, Givens-light e qr e pinv
% di Matlab.
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 4:2:20;

t_sol  = zeros(length(dimensions),5);
e_sol  = zeros(length(dimensions),5);
e_orth = zeros(length(dimensions),3);
e_fact = zeros(length(dimensions),3);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('sottodiagonali: %d\n',n);
    
    dim_m = 150;
    dim_n = 100;
	A = rand(dim_m,dim_n);
    for j = 1:dim_n
        for k = j+n+1:dim_n
            A(k,j) = 0;
        end
    end
    n_cond(i) = cond(A);
    
    sol = ones(dim_n,1);
    norm_sol = norm(sol);
    b = A*sol;

    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,1) = toc;
    e_sol(i,1) = norm(x-sol)/norm_sol;
    e_orth(i,1) = norm(Q*Q'-eye(dim_m));
    e_fact(i,1) = norm(A-Q*R);
    
    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,2) = toc;
    e_sol(i,2) = norm(x-sol)/norm_sol;
    e_orth(i,2) = norm(Q*Q'-eye(dim_m));
    e_fact(i,2) = norm(A-Q*R);
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,3) = toc;
    e_sol(i,3) = norm(x-sol)/norm_sol;
    e_orth(i,3) = norm(Q*Q'-eye(dim_m));
    e_fact(i,3) = norm(A-Q*R);
    
    tic;
    x = pinv(A)*b;
    t_sol(i,4) = toc;
    e_sol(i,4) = norm(x-sol)/norm_sol;
     
    tic;
    x = A\b;
    t_sol(i,5) = toc;
    e_sol(i,5) = norm(x-sol)/norm_sol;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol(:,1), 'b--', ...
         dimensions, e_sol(:,2), 'r--', ...
         dimensions, e_sol(:,3), 'g', ...
         dimensions, e_sol(:,4), 'g:', ...
         dimensions, e_sol(:,5), 'g--');
title("errori di soluzione relativi");
legend('householder-light','givens-light','qr','pinv','backslash');

figure(2); 
semilogy(dimensions, t_sol(:,1), 'b--', ...
         dimensions, t_sol(:,2), 'r--', ...
         dimensions, t_sol(:,3), 'g', ...
         dimensions, t_sol(:,4), 'g:', ...
         dimensions, t_sol(:,5), 'g--');
title("tempi di risoluzione");
legend('householder-light','givens-light','qr','pinv','backslash');

figure(3); 
semilogy(dimensions, e_orth(:,1), 'b--', ...
         dimensions, e_orth(:,2), 'r--', ...
         dimensions, e_orth(:,3), 'g');
title("errori di ortogonalizzazione");
legend('householder-light','givens-light','qr');

figure(4); 
semilogy(dimensions, e_fact(:,1), 'b--', ...
         dimensions, e_fact(:,2), 'r--', ...
         dimensions, e_fact(:,3), 'g');
title("errori di fattorizzazione");
legend('householder-light','givens-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");




% ********************************************************************** %
% ********************************************************************** %
% ***** PARTE 3: sistemi sottodimensionati
% ********************************************************************** %
% ********************************************************************** %

% ********************************************************************** %
% TEST 3.1: : equazioni normali di 2° tipo
% differenza tra Cholesky e pinv con matrici rettangolari orizzontali.
% Sistema visto come A'y = c, con A mxn m>n. Tra le infinite soluzioni, 
% prendo quella di minima norma.
% Sia il sistema A'y = c, la risoluzione è data dal sistema:
% | A'Az = c
% | y = Az
% come prima cosa fattorizzo A'A = R'*R, ottenendo: R'*R*z = c, quindi,
% chiamando t = R*z, ottengo il sistema:
% | R'*t = c  ; | t = R'\c
% | R*z = t     | z = R\t
% dopodiché trovo la soluzione y come y = A*z
% ********************************************************************** %


clear variables;

i = 1;

dimensions = 30:10:110;

t_sol  = zeros(length(dimensions),4);
e_sol  = zeros(length(dimensions),4);
e_fact = zeros(length(dimensions),2);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
    m = 120;
	A = rand(m,n);
    n_cond(i) = cond(A);
    
    sol = ones(m,1)*20;
    norm_sol = norm(sol);
    c = A'*sol;

    tic;
    R = mychol(A'*A); 
    t = R'\c;
    z = R\t;
    y = A*z;
    t_sol(i,1) = toc; 
    e_sol(i,1) = norm(y-sol)/norm_sol; 
    e_fact(i,1) = norm(A'*A-R'*R); 
    
    tic;
    R = chol(A'*A); 
    t = R'\c;
    z = R\t;
    y = A*z;
    t_sol(i,2) = toc; 
    e_sol(i,2) = norm(y-sol)/norm_sol; 
    e_fact(i,2) = norm(A'*A-R'*R); 
    
    tic;
    x = pinv(A')*c;
    t_sol(i,3) = toc;
    e_sol(i,3) = norm(x-sol)/norm_sol;    
    
    tic;
    x = A'\c;
    t_sol(i,4) = toc;
    e_sol(i,4) = norm(x-sol)/norm_sol;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol(:,1), 'm', ...
         dimensions, e_sol(:,2), 'm--', ...
         dimensions, e_sol(:,3), 'g:', ...
         dimensions, e_sol(:,4), 'g--');
title("errori di soluzione relativi");
legend('mychol','chol','pinv','backslash');

figure(2); 
semilogy(dimensions, t_sol(:,1), 'm', ...
         dimensions, t_sol(:,2), 'm--', ...
         dimensions, t_sol(:,3), 'g:', ...
         dimensions, t_sol(:,4), 'g--');
title("tempi di risoluzione");
legend('mychol','chol','pinv','backslash');

figure(4); 
semilogy(dimensions, e_fact(:,1), 'm', ...
         dimensions, e_fact(:,2), 'm--');
title("errori di fattorizzazione");
legend('mychol','chol');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");


% ********************************************************************** %
% TEST 3.2: differenza tra qr e pinv con matrici rettangolari orizzontali. 
% Sistema visto come A'y = c, con A mxn m>n. Tra le infinite soluzioni,  
% prendo quella di minima norma.
% Sia il sistema A'y = c, la risoluzione è data dal sistema:
% | A'Az = c
% | y = Az
% come prima cosa fattorizzo A'A = Q*R, ottenendo: Q*R*z = c, quindi,
% chiamando t = R*z, ottengo il sistema:
% | Q*t = c  ; | t = Q'*c
% | R*z = t    | z = R\t
% dopodiché trovo la soluzione y come y = A*z
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 300:20:500;

t_sol  = zeros(length(dimensions),4);
e_sol  = zeros(length(dimensions),4);
e_fact = zeros(length(dimensions),2);
e_orth = zeros(length(dimensions),2);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
    m = 510;
	A = rand(m,n);
    n_cond(i) = cond(A);
    
    sol = ones(m,1);
    norm_sol = norm(sol);
    c = A'*sol;
    
    tic;
    [Q, R] = myqr(A'*A, 'householder-light'); 
    t = Q'*c;
    z = R\t;
    y = A*z;
    t_sol(i,1) = toc; 
    e_sol(i,1) = norm(y-sol)/norm_sol; 
    e_fact(i,1) = norm(A'*A-Q*R); 
    e_orth(i,1) = norm(Q*Q'-eye(n)); 
    
    tic;
    [Q, R] = qr(A'*A); 
    t = Q'*c;
    z = R\t;
    y = A*z;
    t_sol(i,2) = toc; 
    e_sol(i,2) = norm(y-sol)/norm_sol; 
    e_fact(i,2) = norm(A'*A-Q*R); 
    e_orth(i,2) = norm(Q*Q'-eye(n)); 
        
    tic;
    x = pinv(A')*c;
    t_sol(i,3) = toc;
    e_sol(i,3) = norm(x-sol)/norm_sol;
        
    tic;
    x = A'\c;
    t_sol(i,4) = toc;
    e_sol(i,4) = norm(x-sol)/norm_sol;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol(:,1), 'b--', ...
         dimensions, e_sol(:,2), 'g', ...
         dimensions, e_sol(:,3), 'g:', ...
         dimensions, e_sol(:,4), 'g--');
title("errori di soluzione relativi");
legend('householder-light','qr','pinv','backslash');

figure(2); 
semilogy(dimensions, t_sol(:,1), 'b--', ...
         dimensions, t_sol(:,2), 'g', ...
         dimensions, t_sol(:,3), 'g:', ...
         dimensions, t_sol(:,4), 'g--');
title("tempi di risoluzione");
legend('householder-light','qr','pinv','backslash');

figure(3); 
semilogy(dimensions, e_orth(:,1), 'b--', ...
         dimensions, e_orth(:,2), 'g');
title("errori di ortogonalizzazione");
legend('householder-light','qr');

figure(4); 
semilogy(dimensions, e_fact(:,1), 'b--', ...
         dimensions, e_fact(:,2), 'g');
title("errori di fattorizzazione");
legend('householder-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");



% ********************************************************************** %
% ********************************************************************** %
% ***** PARTE 4: sistemi malcondizionati
% ********************************************************************** %
% ********************************************************************** %

% ********************************************************************** %
% TEST 4.1 sistema sovradimensionato
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 3:1:9;

t_sol  = zeros(length(dimensions),6);
e_sol  = zeros(length(dimensions),6);
e_fact = zeros(length(dimensions),4);
e_orth = zeros(length(dimensions),2);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
    m=10;
	A_start = hilb(m); % hilb crea una matrice quadrata
    A = A_start(:,1:n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    norm_sol = norm(sol);
    b = A*sol;

    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,1) = toc;
    e_sol(i,1) = norm(x-sol)/norm_sol;
    e_orth(i,1) = norm(Q*Q'-eye(m));
    e_fact(i,1) = norm(A-Q*R);
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,2) = toc;
    e_sol(i,2) = norm(x-sol)/norm_sol;
    e_orth(i,2) = norm(Q*Q'-eye(m));
    e_fact(i,2) = norm(A-Q*R);
    

    tic;
    R = mychol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    t_sol(i,3) = toc;
    e_sol(i,3) = norm(x-sol)/norm_sol;
    e_fact(i,3) = norm(A'*A-R'*R);
    
    tic;
    R = chol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    t_sol(i,4) = toc;
    e_sol(i,4) = norm(x-sol)/norm_sol;
    e_fact(i,4) = norm(A'*A-R'*R);
    
    tic;
    x = pinv(A)*b;
    t_sol(i,5) = toc;
    e_sol(i,5) = norm(x-sol)/norm_sol;
    
    tic;
    x = A\b;
    t_sol(i,6) = toc;
    e_sol(i,6) = norm(x-sol)/norm_sol;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol(:,1), 'b--', ...
         dimensions, e_sol(:,2), 'g', ...
         dimensions, e_sol(:,3), 'm', ...
         dimensions, e_sol(:,4), 'm--', ...
         dimensions, e_sol(:,5), 'g:', ...
         dimensions, e_sol(:,6), 'g--');
title("errori di soluzione relativi");
legend('householder-light','qr','mychol','chol','pinv','backslash');

figure(2); 
semilogy(dimensions, t_sol(:,1), 'b--', ...
         dimensions, t_sol(:,2), 'g', ...
         dimensions, t_sol(:,3), 'm', ...
         dimensions, t_sol(:,4), 'm--', ...
         dimensions, t_sol(:,5), 'g:', ...
         dimensions, t_sol(:,6), 'g--');
title("tempi di risoluzione");
legend('householder-light','qr','mychol','chol','pinv','backslash');

figure(3); 
semilogy(dimensions, e_orth(:,1), 'b--', ...
         dimensions, e_orth(:,2), 'g');
title("errori di ortogonalizzazione");
legend('householder-light','qr');

figure(4); 
semilogy(dimensions, e_fact(:,1), 'b--', ...
         dimensions, e_fact(:,2), 'g', ...
         dimensions, e_fact(:,3), 'm', ...
         dimensions, e_fact(:,4), 'm--');
title("errori di fattorizzazione");
legend('householder-light','qr','mychol','chol');

figure(5); 
semilogy(dimensions, n_cond); title("numero di condizionamento");


% ********************************************************************** %
% ********************************************************************** %
% ***** PARTE 5: matrici sparse
% ********************************************************************** %
% ********************************************************************** %

% ************************************************************************
% TEST 5.1: analisi algoritmi con matrici sparse
% ************************************************************************

% NON VA: Warning: Matrix is singular to working precision.
% 
% clear variables;
% 
% i = 1;
% 
% dimensions = 50:50:300;
% 
% t_sol  = zeros(length(dimensions),4);
% e_sol  = zeros(length(dimensions),4);
% e_orth = zeros(length(dimensions),3);
% % e_fact = zeros(length(dimensions),3);
% % n_cond = zeros(length(dimensions),1);
% 
% for n = dimensions
%     fprintf('n: %d\n',n);
%     
% 	A = sprand(n,n,0.01);
% %     n_cond(i) = cond(A);
%     
%     sol = ones(n,1);
%     norm_sol = norm(sol);
%     b = A*sol;
%     
%     tic;
%     [Q,R] = myqr(A, "householder-light");
%     y = Q'*b;
%     x = R\y;
%     t_sol(i,1) = toc;
%     e_sol(i,1) = norm(x-sol)/norm_sol;
%     e_orth(i,1) = norm(Q*Q'-eye(n));
% %     e_fact(i,1) = norm(A-Q*R);
% 
%     tic;
%     [Q,R] = myqr(A, "givens-light");
%     y = Q'*b;
%     x = R\y;
%     t_sol(i,2) = toc;
%     e_sol(i,2) = norm(x-sol)/norm_sol;
%     e_orth(i,2) = norm(Q*Q'-eye(n));
% %     e_fact(i,2) = norm(A-Q*R);
%     
%     tic;
%     [Q,R] = qr(A);
%     y = Q'*b;
%     x = R\y;
%     t_sol(i,3) = toc;
%     e_sol(i,3) = norm(x-sol)/norm_sol;
%     e_orth(i,3) = norm(Q*Q'-eye(n));
% %     e_fact(i,3) = norm(A-Q*R);
%     
%     tic;
%     x = A\b;
%     t_sol(i,4) = toc;
%     e_sol(i,4) = norm(x-sol)/norm_sol;
%     
%     i = i+1;
% 
% end
% 
% figure(1); 
% semilogy(dimensions, e_sol); title("errori di soluzione relativi");
% legend('householder-light','givens-light','qr','backslash');
% 
% figure(2); 
% semilogy(dimensions, t_sol); title("tempi di risoluzione");
% legend('householder-light','givens-light','qr','backslash');
% 
% figure(3); 
% semilogy(dimensions, e_orth); title("errori di ortogonalizzazione");
% legend('householder-light','givens-light','qr');
% 
% figure(4); 
% semilogy(dimensions, e_fact); title("errori di fattorizzazione");
% legend('householder-light','givens-light','qr');



% ************************************************************************
% TEST
% ************************************************************************






    




