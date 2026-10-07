% Monte Carlo BER simulation of RIS-assisted received space shift keying.

clear all;
clc;
tic;

%% Parameters
SNR = -10:2:10;          % SNR in dB
SNRL = 10.^(SNR ./ 10);  % Linear SNR
K = length(SNR);         % Number of SNR points
Nt = 1;                  % Transmit antennas
N = 4;                   % RIS reflecting elements
Nr = 16;                 % Receive antennas

%% Bit Mapping
bpcu = log2(Nr);
bits_total = bpcu * 6e4;                % Total transmitted bits
bit_seq = randi(2, 1, bits_total) - 1;  % Random binary data

%% Monte Carlo Simulation
decoded = zeros(1, bits_total);  % Decoded bit sequence

for ee = 1:K
    sigma = sqrt(((1) / (2 * SNRL(ee))));  % Complex noise variance: 2 * sigma^2
    kk = 1;
    for iteration = 1:bits_total / bpcu
        bits = bit_seq(kk:kk + bpcu - 1);

        % Map bits to a receive antenna.
        dec = sum(bits(1:1 + log2(Nr) - 1) .* [2.^(log2(Nr) - 1:-1:0)]) + 1;
        xu = dec;

        % Generate Rayleigh fading channels and noise.
        H = (randn(N, Nt) + 1i * randn(N, Nt)) / sqrt(2);  % Tx-to-RIS channel
        G = (randn(Nr, N) + 1i * randn(Nr, N)) / sqrt(2);  % RIS-to-Rx channel
        n_B = (sigma) .* (randn(Nr, 1) + 1i * randn(Nr, 1));

        % RIS phase adjustment.
        Phi = exp(-1i * angle(G(xu, :) * H));

        % Received signal.
        y_B = G .* Phi * H + n_B;

        % ML detection.
        metrics = zeros(1, Nr);
        aa = 1;

        for count1 = 1:Nr
            Phi = exp(-1i * angle(G(count1, :) * H));
            GGG = G * Phi * H;
            metrics(aa) = norm(y_B - GGG * 1)^2;
            aa = aa + 1;
        end

        [a1, a2] = min(metrics);

        % Recover antenna bits (MSB first).
        dec1 = zeros(1, log2(Nr));
        for pp = 1:log2(Nr)
            dec1(pp) = mod(floor((a2 - 1) / (2^(log2(Nr) - pp))), 2);
        end

        decoded(kk:kk + log2(Nr) - 1) = dec1;
        kk = kk + log2(Nr);
    end

    % Estimate BER.
    bit_errors = sum(xor(bit_seq, decoded));
    BER(ee) = bit_errors / bits_total;
    disp(BER)
end

%% BER Results
toc;
semilogy(SNR, BER, '-bp', 'LineWidth', 1);
grid on;
xlabel('SNR-dB');
ylabel('BER');
hold on;
