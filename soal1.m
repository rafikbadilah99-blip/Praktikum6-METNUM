% ===========================
% RAFIK BADILAH_L0325009
% PRAKTIKUM-6 METODE NUMERIK
% ===========================

% =========================================================================
% SOAL 1: MENCARI HIMPUNAN PENYELESAIAN SPL
% Metode: Gauss, Gauss-Jordan, LU, Jacobi, dan Gauss-Seidel
% =========================================================================


% Matriks Koefisien A dan Vektor b
A = [ 5 -1 -1  0  0;
     -1  5 -1 -1  0;
     -1 -1  4 -1 -1;
      0  0  1  4 -2;
      0  1 -1  1  4];

b = [-1; 2; 6; 2; -1];
n = length(b);

tol = 0.001;
max_iter = 20;

fprintf('=======================================================\n');
fprintf('                PENYELESAIAN SOAL 1                    \n');
fprintf('=======================================================\n\n');

%% 1A. METODE ELIMINASI GAUSS
Aug = [A b];
for i = 1:n-1
    for j = i+1:n
        factor = Aug(j,i) / Aug(i,i);
        Aug(j,i:n+1) = Aug(j,i:n+1) - factor * Aug(i,i:n+1);
    end
end
x_gauss = zeros(n,1);
x_gauss(n) = Aug(n,n+1) / Aug(n,n);
for i = n-1:-1:1
    x_gauss(i) = (Aug(i,n+1) - Aug(i,i+1:n)*x_gauss(i+1:n)) / Aug(i,i);
end

fprintf('1. Hasil Metode Gauss (p, q, r, s, t):\n');
disp(x_gauss');

%% 1B. METODE GAUSS-JORDAN
AugGJ = [A b];
for i = 1:n
    AugGJ(i,:) = AugGJ(i,:) / AugGJ(i,i);
    for j = 1:n
        if j ~= i
            AugGJ(j,:) = AugGJ(j,:) - AugGJ(j,i) * AugGJ(i,:);
        end
    end
end
x_gj = AugGJ(:, end);

fprintf('2. Hasil Metode Gauss-Jordan (p, q, r, s, t):\n');
disp(x_gj');

%% 1C. DEKOMPOSISI LU
[L, U] = lu(A);
y = L \ b;
x_lu = U \ y;

fprintf('3. Hasil Dekomposisi LU (p, q, r, s, t):\n');
disp(x_lu');

%% 1D. METODE ITERASI JACOBI
C_j = -A;
d_j = zeros(n, 1);
for i = 1:n
    C_j(i,i) = 0;
    C_j(i,:) = C_j(i,:) / A(i,i);
    d_j(i,1) = b(i) / A(i,i);
end
x_jacobi = zeros(n, 1); % Initial guess [0 0 0 0 0]
k_j = 1;
while (k_j <= max_iter)
    xbaru = C_j * x_jacobi + d_j;
    if max(abs(xbaru - x_jacobi)) <= tol
        x_jacobi = xbaru;
        break;
    end
    x_jacobi = xbaru;
    k_j = k_j + 1;
end

fprintf('4. Hasil Iterasi Jacobi (Iterasi ke-%d):\n', k_j);
disp(x_jacobi');

%% 1E. METODE ITERASI GAUSS-SEIDEL
C_gs = -A;
d_gs = zeros(n, 1);
for i = 1:n
    C_gs(i,i) = 0;
    d_gs(i,1) = b(i) / A(i,i);
    C_gs(i,:) = C_gs(i,:) / A(i,i);
end
x_gs = zeros(n, 1); % Initial guess [0 0 0 0 0]
k_gs = 1;
while (k_gs <= max_iter)
    xlama = x_gs;
    for j = 1:n
        x_gs(j) = C_gs(j,:) * x_gs + d_gs(j);
    end
    if max(abs(xlama - x_gs)) <= tol
        break;
    end
    k_gs = k_gs + 1;
end

fprintf('5. Hasil Iterasi Gauss-Seidel (Iterasi ke-%d):\n', k_gs);
disp(x_gs');
