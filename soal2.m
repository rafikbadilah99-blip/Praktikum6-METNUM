% ===========================
% RAFIK BADILAH_L0325009
% PRAKTIKUM-6 METODE NUMERIK
% ===========================

% =========================================================================
% SOAL 2: MENENTUKAN GALAT ABSTRAK (abs(nilai eksak - nilai hampiran))
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

% Nilai Eksak dari Metode Analitik (Hasil Soal 1)
x_eksak = [0; 1; 2; 1; -1];

%% 1. HITUNG ITERASI JACOBI
C_j = -A; d_j = zeros(n, 1);
for i = 1:n
    C_j(i,i) = 0; C_j(i,:) = C_j(i,:) / A(i,i); d_j(i,1) = b(i) / A(i,i);
end
x_jacobi = zeros(n, 1);
k_j = 1;
while (k_j <= max_iter)
    xbaru = C_j * x_jacobi + d_j;
    if max(abs(xbaru - x_jacobi)) <= tol
        x_jacobi = xbaru; break;
    end
    x_jacobi = xbaru; k_j = k_j + 1;
end

%% 2. HITUNG ITERASI GAUSS-SEIDEL
C_gs = -A; d_gs = zeros(n, 1);
for i = 1:n
    C_gs(i,i) = 0; d_gs(i,1) = b(i) / A(i,i); C_gs(i,:) = C_gs(i,:) / A(i,i);
end
x_gs = zeros(n, 1);
k_gs = 1;
while (k_gs <= max_iter)
    xlama = x_gs;
    for j = 1:n, x_gs(j) = C_gs(j,:) * x_gs + d_gs(j); end
    if max(abs(xlama - x_gs)) <= tol, break; end
    k_gs = k_gs + 1;
end

%% 3. PERHITUNGAN GALAT ABSTRAK (abs(nilai eksak - nilai hampiran))
galat_jacobi = abs(x_eksak - x_jacobi);
galat_gs = abs(x_eksak - x_gs);

fprintf('=======================================================\n');
fprintf('                PENYELESAIAN SOAL 2                    \n');
fprintf('=======================================================\n\n');

fprintf('A. Galat Abstrak Metode Jacobi (Iterasi ke-%d):\n', k_j);
fprintf('   - Galat p : %.6f\n', galat_jacobi(1));
fprintf('   - Galat q : %.6f\n', galat_jacobi(2));
fprintf('   - Galat r : %.6f\n', galat_jacobi(3));
fprintf('   - Galat s : %.6f\n', galat_jacobi(4));
fprintf('   - Galat t : %.6f\n\n', galat_jacobi(5));

fprintf('B. Galat Abstrak Metode Gauss-Seidel (Iterasi ke-%d):\n', k_gs);
fprintf('   - Galat p : %.6f\n', galat_gs(1));
fprintf('   - Galat q : %.6f\n', galat_gs(2));
fprintf('   - Galat r : %.6f\n', galat_gs(3));
fprintf('   - Galat s : %.6f\n', galat_gs(4));
fprintf('   - Galat t : %.6f\n', galat_gs(5));
