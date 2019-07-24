% TESINA - Computational Mathematics, Carlo Cabras - 60/73/65113

% da testare:
%  - al variare di m ed n
%     - per ogni algoritmo mio e confronto con MATLAB
%        - tempi di calcolo
%        - condizionamento
%        - norma della differenza con la soluzione vera
% usare pseudoinversa
% 

% IMPORTANTE: lanciare i test uno alla volta per evitare blocchi

% ********************************************************************** %
% ********************************************************************** %
% ***** PARTE 1: matrici quadrate
% ********************************************************************** %
% ********************************************************************** %

% ********************************************************************** %
% TEST 1.1: confronto tra tutti gli algoritmi. A nxn, n max 100
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 50:10:100;

t_sol  = zeros(length(dimensions),7);
e_sol  = zeros(length(dimensions),7);
e_fact = zeros(length(dimensions),7);
e_orth = zeros(length(dimensions),5);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(n,n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    tic;
    [Q,R] = myqr(A, "householder");
    y = Q'*b;
    x = R\y;
    t_sol(i,k)  = toc;
    e_sol(i,k)  = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(n));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(n));
    e_fact(i,k) = norm(A-Q*R);
    k = k+1;
    
    tic;
    [Q,R] = myqr(A, "givens");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(n));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(n));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(n));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;

        
    % Cholesky con n >100 non funziona più, mi dice che la matrice deve
    % essere definita positiva. Probabilmente è dovuto al fatto che A'*A è
    % una matrice di ordine 10'000 e ci sono troppi errori
    
    tic;
    R = mychol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_fact(i,k) = norm(A'*A-R'*R);
    k = k + 1;

    tic;
    R = chol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_fact(i,k) = norm(A'*A-R'*R);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol); title("errori di soluzione");
legend('householder','householder-light','givens','givens-light','qr','mychol','chol');

figure(2); 
semilogy(dimensions, t_sol); title("tempi di risoluzione");
legend('householder','householder-light','givens','givens-light','qr','mychol','chol');

figure(3); 
semilogy(dimensions, e_orth); title("errori di ortogonalizzazione");
legend('householder','householder-light','givens','givens-light','qr');

figure(4); 
semilogy(dimensions, e_sol); title("errori di fattorizzazione");
legend('householder','householder-light','givens','givens-light','qr','mychol','chol');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");


% ********************************************************************** %
% TEST 1.2: visto che Cholesky non funziona più per n >120, in questo
% test non lo uso. Inoltre, visti i tempi più alti per Householder e
% Givens, uso solo le loro versioni light. Così posso usare n più grandi.
% A matrice quadrata nxn
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 50:50:500;

t_sol  = zeros(length(dimensions),3);
e_sol  = zeros(length(dimensions),3);
e_orth = zeros(length(dimensions),3);
e_fact = zeros(length(dimensions),3);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(n,n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(n));
    e_fact(i,k) = norm(A-Q*R);
    k = k+1;

    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(n));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(n));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol); title("errori di soluzione");
legend('householder-light','givens-light','qr');

figure(2); 
semilogy(dimensions, t_sol); title("tempi di risoluzione");
legend('householder-light','givens-light','qr');

figure(3); 
semilogy(dimensions, e_orth); title("errori di ortogonalizzazione");
legend('householder-light','givens-light','qr');

figure(4); 
semilogy(dimensions, e_sol); title("errori di fattorizzazione");
legend('householder-light','givens-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");


% ********************************************************************** %
% TEST 1.3: comparazione tra mio Householder-light e qr di MATLAB, con
% n grandi. A matrice quadrata nxn
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 1000:500:1500;

t_sol  = zeros(length(dimensions),2);
e_sol  = zeros(length(dimensions),2);
e_orth = zeros(length(dimensions),2);
e_fact = zeros(length(dimensions),2);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(n,n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(n));
    e_fact(i,k) = norm(A-Q*R);
    k = k+1;

    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(n));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol); title("errori di soluzione");
legend('householder-light','qr');

figure(2); 
semilogy(dimensions, t_sol); title("tempi di risoluzione");
legend('householder-light','qr');

figure(3); 
semilogy(dimensions, e_orth); title("errori di ortogonalizzazione");
legend('householder-light','qr');

figure(4); 
semilogy(dimensions, e_sol); title("errori di fattorizzazione");
legend('householder-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");

% ********************************************************************** %
% TEST 1.4: matrice triangolare superiore con diagonali in più
% dimensioni di A fissate a 500x500, cambia il numero delle sottodiagonali.
% Confrontro tra Householder-light, Givens-light e qr di MATLAB.
% ********************************************************************** %
% funziona male per matrici quasi triangolari, perché mi dà l'errore:
% "Warning: Matrix is close to singular or badly scaled. Results may be 
% inaccurate. RCOND =  6.550039e-22." quindi parto con 4 sottodiagonali.
% Si nota comunque che Givens-light vince su Householder-light, ma non su
% qr di Matlab

clear variables;

i = 1;

dimensions = 4:2:20;

t_sol  = zeros(length(dimensions),3);
e_sol  = zeros(length(dimensions),3);
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
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(dim));
    e_fact(i,k) = norm(A-Q*R);
    k = k+1;
    
    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(dim));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(dim));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol); title("errori di soluzione");
legend('householder-light','givens-light','qr');

figure(2); 
semilogy(dimensions, t_sol); title("tempi di risoluzione");
legend('householder-light','givens-light','qr');

figure(3); 
semilogy(dimensions, e_orth); title("errori di ortogonalizzazione");
legend('householder-light','givens-light','qr');

figure(4); 
semilogy(dimensions, e_sol); title("errori di fattorizzazione");
legend('householder-light','givens-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");


% ********************************************************************** %
% ********************************************************************** %
% ***** PARTE 2: matrici rettangolari verticali
% ********************************************************************** %
% ********************************************************************** %

% ********************************************************************** %
% TEST 2.1: differenza tra Cholesky con matrici rettangolari verticali
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 10:10:100;

t_sol = zeros(length(dimensions),2);
e_sol = zeros(length(dimensions),2);
e_fact = zeros(length(dimensions),2);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(110,n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando

    tic;
    R = mychol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_fact(i,k) = norm(A'*A-R'*R);
    k = k + 1;
    
    tic;
    R = chol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_fact(i,k) = norm(A'*A-R'*R);
    k = k + 1;
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol); title("errori di soluzione");
legend('mychol','chol');

figure(2); 
semilogy(dimensions, t_sol); title("tempi di risoluzione");
legend('mychol','chol');

figure(4); 
semilogy(dimensions, e_sol); title("errori di fattorizzazione");
legend('mychol','chol');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");


% ********************************************************************** %
% TEST 2.2: differenza tra myqr e qr di MATLAB. A matrice rettangolare 
% m x n con m > n.
% ********************************************************************** %
% Uso la mia versione di qr con householder-light, comparandola alla qr di
% Matlab

clear variables;

i = 1;

dimensions = 1000:500:1500;

t_sol  = zeros(length(dimensions),2);
e_sol  = zeros(length(dimensions),2);
e_orth = zeros(length(dimensions),2);
e_fact = zeros(length(dimensions),2);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
    m = 1600;
	A = rand(m,n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(m));
    e_fact(i,k) = norm(A-Q*R);
    k = k+1;

    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(m));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol); title("errori di soluzione");
legend('householder-light','qr');

figure(2); 
semilogy(dimensions, t_sol); title("tempi di risoluzione");
legend('householder-light','qr');

figure(3); 
semilogy(dimensions, e_orth); title("errori di ortogonalizzazione");
legend('householder-light','qr');

figure(4); 
semilogy(dimensions, e_sol); title("errori di fattorizzazione");
legend('householder-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");

% ********************************************************************** %
% TEST 2.3: matrice triangolare superiore con diagonali in più. Le
% dimensioni di A  sono fissate a 800x300, cambia il numero delle 
% sottodiagonali. Confrontro tra Householder-light, Givens-light e qr di
% Matlab.
% ********************************************************************** %

clear variables;

i = 1;

dimensions = 4:2:20;

t_sol  = zeros(length(dimensions),3);
e_sol  = zeros(length(dimensions),3);
e_orth = zeros(length(dimensions),3);
e_fact = zeros(length(dimensions),3);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('sottodiagonali: %d\n',n);
    
    dim_m = 400;
    dim_n = 100;
	A = rand(dim_m,dim_n);
    for j = 1:dim_n
        for k = j+n+1:dim_n
            A(k,j) = 0;
        end
    end
    n_cond(i) = cond(A);
    
    sol = ones(dim_n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(dim_m));
    e_fact(i,k) = norm(A-Q*R);
    k = k+1;
    
    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(dim_m));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_orth(i,k) = norm(Q*Q'-eye(dim_m));
    e_fact(i,k) = norm(A-Q*R);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol); title("errori di soluzione");
legend('householder-light','givens-light','qr');

figure(2); 
semilogy(dimensions, t_sol); title("tempi di risoluzione");
legend('householder-light','givens-light','qr');

figure(3); 
semilogy(dimensions, e_orth); title("errori di ortogonalizzazione");
legend('householder-light','givens-light','qr');

figure(4); 
semilogy(dimensions, e_sol); title("errori di fattorizzazione");
legend('householder-light','givens-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");




% ********************************************************************** %
% ********************************************************************** %
% ***** PARTE 3: matrici rettangolari orizzontali
% ********************************************************************** %
% ********************************************************************** %

% DA RIVEDERE COMPLETAMENTE. Qua si usano le equazioni normali di secondo
% tipo, cosa che io non ho fatto manco per il cazzo.

% ********************************************************************** %
% TEST 3.1: differenza tra Cholesky con matrici rettangolari orizzontali.
% ********************************************************************** %
% NON FUNZIONA. Il problema sta negli autovalori ed a problemi di calcolo:
% risulta che il primo autovalore è uno "zero negativo", ovvero poco poco
% più piccolo di zero, es. -2.23e-16. Quindi risulta che la matrice non è
% definita positiva e chol() non può proseguire con i calcoli.

clear variables;

i = 1;

dimensions = 30:10:100;

t_sol  = zeros(length(dimensions),2);
e_sol  = zeros(length(dimensions),2);
e_fact = zeros(length(dimensions),2);
n_cond = zeros(length(dimensions),1);

for m = dimensions
    fprintf('n: %d\n', m);
    
    n = 110;
	A = rand(m,n);
    n_cond(i) = cond(A);
    
    sol = ones(n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando

    tic;
    R = mychol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_fact(i,k) = norm(A'*A-R'*R);
    k = k + 1;
    
    tic;
    R = chol(A'*A); 
    y = R'\(A'*b);
    x = R\y;
    t_sol(i,k) = toc;
    e_sol(i,k) = norm(x-sol);
    e_fact(i,k) = norm(A'*A-R'*R);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol); title("errori di soluzione");
legend('mychol','chol');

figure(2); 
semilogy(dimensions, t_sol); title("tempi di risoluzione");
legend('mychol','chol');

figure(4); 
semilogy(dimensions, e_sol); title("errori di fattorizzazione");
legend('mychol','chol');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");




% ********************************************************************** %
% TEST 3.2: differenza tra qr con matrici rettangolari orizzontali.
% ********************************************************************** %
% tra le infinite soluzioni, prendo quella di minima norma

clear variables;

i = 1;

dimensions = 200:30:480;

t_sol  = zeros(length(dimensions),2);
e_sol  = zeros(length(dimensions),2);
e_fact = zeros(length(dimensions),2);
e_orth = zeros(length(dimensions),2);
n_cond = zeros(length(dimensions),1);

for n = dimensions
    fprintf('n: %d\n',n);
    
    m = 500;
	A = rand(m,n);
    n_cond(i) = cond(A);
    
    sol = ones(m,1);
    c = A'*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    
    tic;
    [Q, R] = myqr(A'*A, 'householder-light'); 
    yy = Q'*c;
    z = R\yy;
    y = A*z;
    t_sol(i,k) = toc; 
    e_sol(i,k) = norm(y-sol); 
    e_fact(i,k) = norm(A'*A-Q*R); 
    e_orth(i,k) = norm(Q*Q'-eye(n)); 
    k = k + 1;
    
    tic;
    [Q, R] = qr(A'*A); 
    yy = Q'*c;
    z = R\yy;
    y = A*z;
    t_sol(i,k) = toc; 
    e_sol(i,k) = norm(y-sol); 
    e_fact(i,k) = norm(A'*A-Q*R); 
    e_orth(i,k) = norm(Q*Q'-eye(n)); 
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, e_sol); title("errori di soluzione");
legend('householder-light','qr');

figure(2); 
semilogy(dimensions, t_sol); title("tempi di risoluzione");
legend('householder-light','qr');

figure(3); 
semilogy(dimensions, e_orth); title("errori di ortogonalizzazione");
legend('householder-light','qr');

figure(4); 
semilogy(dimensions, e_sol); title("errori di fattorizzazione");
legend('householder-light','qr');

figure(5); 
plot(dimensions, n_cond); title("numero di condizionamento");

% RISULTATI DA STUDIARE. errori di soluzione altissimissimi


% ********************************************************************** %
% TEST
% ********************************************************************** %



% ************************************************************************
% TEST
% ************************************************************************



% ************************************************************************
% TEST
% ************************************************************************






    




