% TESINA - Computational Mathematics, Carlo Cabras - 60/73/65113

% idea: lanciare i primi test con tutte le versioni e poi, una volta
% appurato che le versioni base di householder/givens sono scarse,
% usare solo le versioni light

% da testare:
%  - al variare di m ed n
%     - per ogni algoritmo mio e confronto con MATLAB
%        - tempi di calcolo
%        - condizionamento
%        - norma della differenza con la soluzione vera
% 
% 

% **** lanciare i test un blocco alla volta

% ************************************************************************
% TEST 1: confronto tra tutti gli algoritmi. A nxn, n max 100
% ************************************************************************

i = 1;

dimensions = 50:10:100;

tempi = zeros(length(dimensions),7);
differenze = zeros(length(dimensions),7);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(n,n);
    
    sol = ones(n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    tic;
    [Q,R] = myqr(A, "householder");
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k+1;
    
    tic;
    [Q,R] = myqr(A, "givens");
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;
    
    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;

        
    % Cholesky con n >100 non funziona più, mi dice che la matrice deve
    % essere definita positiva
    
    tic;
    R = mychol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;

    tic;
    R = chol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, differenze); title("differenze");
legend('householder','householder-light','givens','givens-light','qr','mychol','chol');

figure(2); 
semilogy(dimensions, tempi); title("tempi");
legend('householder','householder-light','givens','givens-light','qr','mychol','chol');

% ************************************************************************
% TEST 2: visto che Cholesky non funziona più per n >120, in questo
% test non lo uso. Inoltre, visti i tempi più alti per Householder e
% Givens, uso solo le loro versioni light. Così posso usare n più grandi.
% A matrice quadrata nxn
% ************************************************************************

i = 1;

dimensions = 50:50:500;

tempi = zeros(length(dimensions),3);
differenze = zeros(length(dimensions),3);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(n,n);
    
    sol = ones(n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k+1;

    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, differenze); title("differenze");
legend('householder-light','givens-light','qr');

figure(2); 
semilogy(dimensions, tempi); title("tempi");
legend('householder-light','givens-light','qr');

% ************************************************************************
% TEST 3: comparazione tra mio Householder-light e qr di MATLAB, con
% n grandi. A matrice quadrata nxn
% ************************************************************************


i = 1;

dimensions = 1000:500:1500;

tempi = zeros(length(dimensions),2);
differenze = zeros(length(dimensions),2);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(n,n);
    
    sol = ones(n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k+1;

    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, differenze); title("differenze");
legend('householder-light','qr');

figure(2); 
semilogy(dimensions, tempi); title("tempi");
legend('householder-light','qr');

% ************************************************************************
% TEST 4: differenza tra mychol e chol di MATLAB. A matrice quadrata nxn.
% ************************************************************************

i = 1;

dimensions = 10:10:100;

tempi = zeros(length(dimensions),2);
differenze = zeros(length(dimensions),2);

for n = dimensions
    fprintf('n: %d\n',n);
    
	A = rand(n,n);
    
    sol = ones(n,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    
    tic;
    R = mychol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;

    tic;
    R = chol(A'*A);
    y = R'\(A'*b);
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, differenze); title("differenze");
legend('mychol','chol');

figure(2); 
semilogy(dimensions, tempi); title("tempi");
legend('mychol','chol');

% ************************************************************************
% TEST 5: matrice triangolare superiore con diagonali in più
% dimensioni di A fissate a 500x500, cambia il numero delle sottodiagonali.
% Confrontro tra Householder-light, Givens-light e qr di MATLAB.
% ************************************************************************

i = 1;

dimensions = 2:2:20;

tempi = zeros(length(dimensions),3);
differenze = zeros(length(dimensions),3);

for n = dimensions
    fprintf('sottodiagonali: %d\n',n);
    
    dim = 300;
	A = rand(dim,dim);
    for j = 1:dim
        for k = j+n+1:dim
            A(k,j) = 0;
        end
    end
    
    sol = ones(dim,1);
    b = A*sol;
    
    k = 1; % indice nelle matrici dei risultati; comodo perché non sto a
           % tenere conto di quanti algoritmi sto testando
    
    tic;
    [Q,R] = myqr(A, "householder-light");
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k+1;

    tic;
    [Q,R] = myqr(A, "givens-light");
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;
    
    tic;
    [Q,R] = qr(A);
    y = Q'*b;
    x = R\y;
    tempi(i,k) = toc;
    differenze(i,k) = norm(x-sol);
    k = k + 1;
    
    i = i+1;

end

figure(1); 
semilogy(dimensions, differenze); title("differenze");
legend('householder-light','givens-light','qr');

figure(2); 
semilogy(dimensions, tempi); title("tempi");
legend('householder-light','givens-light','qr');


% ************************************************************************
% TEST 6: problemi duali
% ************************************************************************


% ************************************************************************
% TEST 7:
% ************************************************************************



    




