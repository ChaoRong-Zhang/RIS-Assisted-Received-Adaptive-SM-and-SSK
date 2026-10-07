% Monte Carlo BER simulation of received generalized spatial modulation.

clc;
clear all;
format shortE
tic;

%% Parameters
SNR = -10:2:10;          % SNR in dB
SNRL = 10.^(SNR ./ 10);  % Linear SNR
K = length(SNR);         % Number of SNR points
Nt = 1;                  % Transmit antennas
N = 4;                   % RIS reflecting elements
Nr = 16;                 % Receive antennas
ND = 5;                  % Candidate receive antennas
M = 2;                   % Modulation order
Ns = 2;                  % Selected antennas per combination
Es = Ns;

%% Antenna Combinations
cmb = nchoosek(ND, Ns);
Txindx = nchoosek(1:ND, Ns);
p = floor(log2(cmb));  % Antenna-combination bits
q = log2(M);           % Constellation bits
Nc = 2.^p;             % Usable antenna combinations
bpcu = p + q;

%% Bit Mapping
bits_total = bpcu * 6e4;                % Total transmitted bits
bit_seq = randi(2, 1, bits_total) - 1;  % Random binary data

%% Constellation
ss = qammod(0:M - 1, M, 'gray', 'UnitAveragePower', true);
ss = ss ./ sqrt(mean(abs(ss).^2));
if M == 2
    ss = ss .* exp(1i * pi / 4);  % Rotate BPSK by pi/4
end

%% Monte Carlo Simulation
decoded = zeros(1, bits_total);  % Decoded bit sequence

for ee = 1:K
    sigma = sqrt(((1) / (2 * SNRL(ee))));  % Complex noise variance: 2 * sigma^2
    kk = 1;
    for iteration = 1:bits_total / bpcu
        bits = bit_seq(kk:kk + bpcu - 1);

        % Map bits to an antenna combination and a symbol.
        CbNo = sum(bits(1:1 + p - 1) .* [2.^(p - 1:-1:0)]) + 1;
        SymNo = sum(bits(p + 1:end) .* [2.^(log2(M) - 1:-1:0)]) + 1;
        s = ss(SymNo);

        % Generate Rayleigh fading channels.
        H = (randn(N, Nt) + 1i * randn(N, Nt)) / sqrt(2);  % Tx-to-RIS channel
        G = (randn(Nr, N) + 1i * randn(Nr, N)) / sqrt(2);  % RIS-to-Rx channel

        % Assign RIS elements to the selected antennas.
        h = G(Txindx(CbNo, :), :);
        G2 = zeros(1, N);
        for t = 1:N
            G2(1, t) = h(floor(mod(t, Ns) + 1), t);
        end

        % RIS phase adjustment.
        Phi = exp(-1i * angle(G2 * H));

        % Received signal.
        n_B = (sigma) .* (randn(Nr, 1) + 1i * randn(Nr, 1));
        A1 = G * Phi * H;
        y_B = sqrt(Es) * A1 * s + n_B;

        % ML detection.
        metrics = zeros(1, 2^bpcu);
        index1 = metrics;
        index2 = metrics;
        aa = 1;

        for count1 = 1:Nc
            h_3 = G(Txindx(count1, :), :);
            G3 = zeros(1, N);
            for t = 1:N
                G3(1, t) = h_3(floor(mod(t, Ns) + 1), t);
            end
            Phi1 = exp(-1i * angle(G3 * H));
            A3 = G * Phi1 * H;
            for count2 = 1:M
                metrics(aa) = norm((y_B - sqrt(Es) * A3 * ss(count2)), 2);
                index1(aa) = count1;
                index2(aa) = count2;
                aa = aa + 1;
            end
        end

        [a5, a6] = min(metrics);
        a2 = index1(a6);
        a4 = index2(a6);

        % Recover antenna-combination and symbol bits (MSB first).
        dec_IM_Rx = zeros(1, p);
        for pp = 1:p
            dec_IM_Rx(pp) = mod(floor((a2 - 1) / (2^(p - pp))), 2);
        end

        dec_symbol_Rx = zeros(1, log2(M));
        for pp = 1:log2(M)
            dec_symbol_Rx(pp) = mod(floor((a4 - 1) / (2^(log2(M) - pp))), 2);
        end

        decoded(kk:kk + p + log2(M) - 1) = [dec_IM_Rx, dec_symbol_Rx];
        kk = kk + p + log2(M);
    end

    % Estimate BER.
    bit_errors = sum(xor(bit_seq, decoded));
    BER(ee) = bit_errors / bits_total;
    disp(BER)
end

%% BER Results
toc;
semilogy(SNR, BER, '-mp', 'LineWidth', 1);
grid on;
xlabel('SNR-dB');
ylabel('BER');
hold on;
