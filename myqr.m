function [Q, A] = myqr(A, triang_type)

[m, n] = size(A);
Q = eye(m);

if triang_type == "householder" % testato e funzionante
    
    if m > n
        nf = n;
    else
        nf = m-1;
    end
    
    for j = 1:nf
        x = A(j:end, j);
        [~, ~, Ht] = housemat(x); % Ht = H tilde
        H = [ eye(j-1)            zeros(j-1, m-j+1); ...
              zeros(m-j+1, j-1)   Ht                ];
        A = H*A;
        Q = Q*H;
    end

elseif triang_type == "householder-light" % testato e funzionante
    % in questa versione evito di calcolare la matrice H: così le
    % moltiplicazioni avvengono solo tra matrice e vettore
    
    if m > n
        nf = n;
    else
        nf = m-1;
    end
    
    for j = 1:nf
        x = A(j:end, j);
        % [~, ~, Ht] = housemat(x); % Ht = H tilde
        [w, ~, ~] = housemat(x);

        % H = [ eye(j-1)            zeros(j-1, m-j+1); ...
        %       zeros(m-j+1, j-1)   Ht                ];
        
        % A = H*A;
        A = H_per_A(A, w, j, m-j+1);
        % Q = Q*H;
        Q = Q_per_H(Q, w, j);
    end
  
    
elseif triang_type == "givens" % testato e funzionante

    if m > n
        nf = n;
    else
        nf = m-1;
    end
    
    for j = 1:nf
        for i = j+1:m
            x = A(1:end, j); % dentro il ciclo perché cambia ogni volta!
            
            [~, ~, G] = mygivens(x, j, i);
            
            A = G*A;
            Q = Q*G';
        end
    end
        
elseif triang_type == "givens-light"
    C = zeros(m,n);
    S = zeros(m,n);
    
    if m > n
        nf = n;
    else
        nf = m-1;
    end
    
    for k = 1:nf
        for i = k+1:m
%             [i j]
            x = A(1:end, k); % dentro il ciclo perché cambia ogni volta!
            if x(i) ~= 0
                [c, s, ~] = mygivens(x, k, i);
                for j = k:n
    %             A = G*A;
                    t = A(k,j)*c + s*A(i,j);
                    A(i,j) = -s*A(k,j) + c*A(i,j);
                    A(k,j) = t;


                end           
                % ricordando che k: numero da modificare, i: numero da 
                % azzerare;
                % per fare Q = Q*G', si devono modificare le colonne k e i
                % della matrice Q, qunidi uso l'indice j per scorrere tra 
                % le righe
                for j = 1:m
                    t = Q(j,k)*c + s*Q(j,i);
                    Q(j,i) = -s*Q(j,k) + c*Q(j,i);
                    Q(j,k) = t;

                end
            end
        end
    end
    
else
    error("tipo sbagliato, dev'essere 'householder','householder-light', 'givens', 'givens-light");
end

