function compare_norm_generalH_FIXED()

% =========================
% Sizes
% =========================
Nr = 3;
P  = 3;
Q  = 3;
n  = 4;

% =========================
% Random complex data
% =========================
beta = randn(Nr,1) + 1i*randn(Nr,1);
h    = randn(Nr,P) + 1i*randn(Nr,P);
g    = randn(P,Q)  + 1i*randn(P,Q);

% =========================
% General matrices H_{p,q}
% =========================
H = cell(P,Q);
for p = 1:P
    for q = 1:Q
        H{p,q} = randn(n) + 1i*randn(n);
    end
end

% =========================
% 1️⃣ Direct A
% =========================
A = zeros(n);

for r = 1:Nr
    for p = 1:P
        for q = 1:Q
            A = A + beta(r)*h(r,p)*h(r,q)*g(p,q)*H{p,q};
        end
    end
end

norm_direct = norm(A,'fro')^2;

% =========================
% 2️⃣ Fully expanded (CORRECT conjugation)
% =========================
norm_expanded = 0;

for r = 1:Nr
    for s = 1:Nr
        for p = 1:P
            for q = 1:Q
                for u = 1:P
                    for v = 1:Q

                        inner = trace(H{p,q}' * H{u,v});

                        term = conj(beta(r))*beta(s) ...
                             * conj(h(r,p)*h(r,q)*g(p,q)) ...
                             * (h(s,u)*h(s,v)*g(u,v)) ...
                             * inner;

                        norm_expanded = norm_expanded + term;

                    end
                end
            end
        end
    end
end

norm_expanded = real(norm_expanded);

% =========================
% 3️⃣ Vectorized with Gram matrix
% =========================
K = P*Q;
index = @(p,q) (p-1)*Q + q;

Gram = zeros(K,K);
for p = 1:P
    for q = 1:Q
        k = index(p,q);
        for u = 1:P
            for v = 1:Q
                l = index(u,v);
                Gram(k,l) = trace(H{p,q}' * H{u,v});
            end
        end
    end
end

M = zeros(K,Nr);
for r = 1:Nr
    for p = 1:P
        for q = 1:Q
            k = index(p,q);
            M(k,r) = h(r,p)*h(r,q)*g(p,q);
        end
    end
end

norm_vectorized = 0;

for r = 1:Nr
    for s = 1:Nr
        norm_vectorized = norm_vectorized ...
            + conj(beta(r))*beta(s) ...
            * (M(:,r)' * Gram * M(:,s));
    end
end

norm_vectorized = real(norm_vectorized);

% =========================
% Results
% =========================
fprintf('Direct        : %.12f\n', norm_direct);
fprintf('Fully expanded: %.12f\n', norm_expanded);
fprintf('Vectorized    : %.12f\n', norm_vectorized);

fprintf('\nDifferences:\n');
fprintf('|direct - expanded|   = %.3e\n', abs(norm_direct - norm_expanded));
fprintf('|direct - vectorized| = %.3e\n', abs(norm_direct - norm_vectorized));

end
