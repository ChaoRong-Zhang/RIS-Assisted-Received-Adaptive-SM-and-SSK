% Spectral-efficiency comparison of RIS-assisted index modulation schemes.

clear all;
clc

%% Parameters
Nt = 1;                      % Transmit antennas
N = 8;                       % RIS reflecting elements
Nr = 16;                     % Receive antennas
Nd = [2, 4, 8, 16, 32, 64];  % Candidate receive-antenna counts
M = 2;                       % Modulation order
Ns = 2;                      % Selected antennas per combination

%% Result Arrays
bpcu_RASM = zeros(1, length(Nd));
bpcu_RSM = zeros(1, length(Nd));
bpcu_RGSM = zeros(1, length(Nd));
bpcu_RGSSK = zeros(1, length(Nd));
bpcu_RASSK = zeros(1, length(Nd));

%% Spectral Efficiency
for ii = 1:length(Nd)
    i = Nd(ii);

    % RASM
    b11 = 2^i - 1;
    b1_RASM = floor(log2(b11));
    b2_RASM = log2(M);
    bpcu_RASM(ii) = b1_RASM + b2_RASM;

    % RSM
    b1_RSM = log2(i);
    b2_RSM = log2(M);
    bpcu_RSM(ii) = b1_RSM + b2_RSM;

    % RGSM
    b12 = nchoosek(i, Ns);
    b1_RGSM = floor(log2(b12));
    b2 = log2(M);
    bpcu_RGSM(ii) = b1_RGSM + b2;

    % RGSSK
    b13 = nchoosek(Nd, Ns);
    bpcu_RGSSK(ii) = floor(log2(b12));

    % RASSK
    b14 = 2^i - 1;
    bpcu_RASSK(ii) = floor(log2(b14));
end

%% Results
disp(bpcu_RASM);
disp(bpcu_RASSK);
disp(bpcu_RSM);
disp(bpcu_RGSM);
disp(bpcu_RGSSK);

%% Plot
figure;
x = Nd;
y1 = bpcu_RASM;
y2 = bpcu_RASSK;
y3 = bpcu_RSM;
y4 = bpcu_RGSM;
y5 = bpcu_RGSSK;
semilogy(x, y1, 'x-r', 'LineWidth', 1.5);
grid on;
hold on;
semilogy(x, y2, 's-b', 'LineWidth', 1.5);
grid on;
hold on;
semilogy(x, y3, 'o--k', 'LineWidth', 1.5);
grid on;
hold on;
semilogy(x, y4, '*-.g', 'LineWidth', 1.5);
grid on;
hold on;
semilogy(x, y5, '>--m', 'LineWidth', 1.5);
grid on;
hold on;
legend('RSAM', 'RASSK', 'RAM', 'RGAM', 'RGSSK', 'location', 'best');
