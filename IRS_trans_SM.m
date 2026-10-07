% Monte Carlo BER simulation of RIS-assisted received spatial modulation.

clear all;
clc;
tic;

%% Parameters
SNR = -10:2:20;          % SNR in dB
SNRL = 10.^(SNR ./ 10);  % Linear SNR
K = length(SNR);         % Number of SNR points
Nt = 1;                  % Transmit antennas
N = 4;                   % RIS reflecting elements
ND = 8;                  % Candidate receive antennas
Nr = 8;                  % Receive antennas
M = 2;                   % Modulation order

%% Bit Mapping
bpcu = log2(Nr) + log2(M);
bits_total = bpcu * 6e4;                % Total transmitted bits
bit_seq = randi(2, 1, bits_total) - 1;  % Random binary data

%% Constellation
ss = qammod(0:M - 1, M, 'Gray');
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

        % Map bits to an antenna and a symbol.
        dec_IM = sum(bits(1:1 + log2(Nr) - 1) .* [2.^(log2(Nr) - 1:-1:0)]) + 1;
        xu = dec_IM;
        dec_symbol = sum(bits(log2(Nr) + 1:end) .* [2.^(log2(M) - 1:-1:0)]) + 1;
        symbol = ss(dec_symbol);

        % Generate Rayleigh fading channels and noise.
        H = (randn(N, Nt) + 1i * randn(N, Nt)) / sqrt(2);  % Tx-to-RIS channel
        G = (randn(ND, N) + 1i * randn(ND, N)) / sqrt(2);  % RIS-to-Rx channel
        n_B = (sigma) .* (randn(ND, 1) + 1i * randn(ND, 1));

        % RIS phase adjustment.
        Phi = exp(-1i * angle(G(xu, :) * H));

        % Received signal.
        HHH = G * Phi * H;
        y_B = HHH * symbol + n_B;

        % ML detection.
        metrics = zeros(1, 2^bpcu);
        index1 = metrics;
        index2 = metrics;
        aa = 1;

        for count1 = 1:Nr
            Phi = exp(-1i * angle(G(count1, :) * H));
            GGG = G * Phi * H;
            for count2 = 1:M
                metrics(aa) = norm(y_B - GGG * ss(count2))^2;
                index1(aa) = count1;
                index2(aa) = count2;
                aa = aa + 1;
            end
        end

        [a5, a6] = min(metrics);
        a2 = index1(a6);
        a4 = index2(a6);

        % Recover antenna and symbol bits (MSB first).
        dec_IM_Rx = zeros(1, log2(Nr));
        for pp = 1:log2(Nr)
            dec_IM_Rx(pp) = mod(floor((a2 - 1) / (2^(log2(Nr) - pp))), 2);
        end

        dec_symbol_Rx = zeros(1, log2(M));
        for pp = 1:log2(M)
            dec_symbol_Rx(pp) = mod(floor((a4 - 1) / (2^(log2(M) - pp))), 2);
        end

        decoded(kk:kk + log2(Nr) + log2(M) - 1) = [dec_IM_Rx, dec_symbol_Rx];
        kk = kk + log2(Nr) + log2(M);
    end

    % Estimate BER.
    bit_errors = sum(xor(bit_seq, decoded));
    BER(ee) = bit_errors / bits_total;
    disp(BER)
end

%% BER Results
toc;
semilogy(SNR, BER, '-yp', 'LineWidth', 1);
grid on;
xlabel('SNR-dB');
ylabel('BER');
hold on;
