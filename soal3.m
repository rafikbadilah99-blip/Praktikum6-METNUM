% ===========================
% RAFIK BADILAH_L0325009
% PRAKTIKUM-6 METODE NUMERIK
% ===========================

% =========================================================================
% SOAL 3: PENGUJIAN RUNNING TIME PROSES (Hanya Proses)
% =========================================================================


% Input Data (Tidak dihitung dalam running time)
A = [ 5 -1 -1  0  0;
     -1  5 -1 -1  0;
     -1 -1  4 -1 -1;
      0  0  1  4 -2;
      0  1 -1  1  4];

b = [-1; 2; 6; 2; -1];
n = length(b);
tol = 0.001;
max_iter = 20;

%% 1. RUNNING TIME - GAUSS
tic; % Mulai ukur algoritma saja
Aug = [A b];
for i = 1:n-1
    for j = i+1:n
        factor = Aug(j,i) / Aug(i,i);
        Aug(j,i:n+1) = Aug(j,i:n+1) - factor * Aug(i,i:n+1);
    end
end
x_g = zeros(n,1);
x_g(n) = Aug(n,n+1) / Aug(n,n);
for i = n-1:-1:1
    x_g(i) = (Aug(i,n+1) - Aug(i,i+1:n)*x_g(i+1:n)) / Aug(i,i);
end
t_gauss = toc; % Selesai ukur

%% 2. RUNNING TIME - GAUSS JORDAN
tic;
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
t_gj = toc;

%% 3. RUNNING TIME - LU
tic;
[L, U] = lu(A);
y = L \ b;
x_lu = U \ y;
t_lu = toc;

%% 4. RUNNING TIME - JACOBI
tic;
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
t_jacobi = toc;

%% 5. RUNNING TIME - GAUSS-SEIDEL
tic;
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
t_gs = toc;

%% TAMPILKAN HASIL WAKTU EKSEKUSI
fprintf('=======================================================\n');
fprintf('                PENYELESAIAN SOAL 3                    \n');
fprintf('=======================================================\n\n');
fprintf('Running Time Algoritma (diukur menggunakan tic-toc):\n');
fprintf('1. Metode Gauss        : %.6f detik\n', t_gauss);
fprintf('2. Metode Gauss-Jordan : %.6f detik\n', t_gj);
fprintf('3. Dekomposisi LU      : %.6f detik\n', t_lu);
fprintf('4. Iterasi Jacobi      : %.6f detik\n', t_jacobi);
fprintf('5. Iterasi Gauss-Seidel: %.6f detik\n', t_gs);
