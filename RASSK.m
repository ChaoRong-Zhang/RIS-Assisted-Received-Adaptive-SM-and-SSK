% Monte Carlo BER simulation of received adaptive space shift keying.

clear all;
clc;
format shortE
tic;

%% Parameters
SNR = -10:2:20;          % SNR in dB
SNRL = 10.^(SNR ./ 10);  % Linear SNR
K = length(SNR);         % Number of SNR points
Nt = 1;                  % Transmit antennas
N = 8;                   % RIS reflecting elements
Nr = 5;                  % Receive antennas
ND = 5;                  % Candidate receive antennas
Ns = 2^ND - 1;           % Nonempty antenna combinations
Nw = sort(randperm(ND));
cmb = 0;
aTx = [];

%% Antenna Combinations
for Nf = 1:ND
    cmbb = nchoosek(ND, Nf);
    cmb = cmb + cmbb;
end

for iii = 1:ND
    Tx1 = nchoosek(1:ND, iii);
    Tx2 = padarray(Tx1, [0, ND - iii], 'post');
    aTx = [aTx; Tx2];
end

%% Bit Mapping
p1 = floor(log2(Ns));  % Antenna-combination bits
Nc = 2.^p1;            % Usable antenna combinations
bpcu = p1;
bits_total = bpcu * 6e4;  % Total transmitted bits

%% Monte Carlo Simulation
bit_seq = randi(2, 1, bits_total) - 1;  % Random binary data
decoded = zeros(1, bits_total);         % Decoded bit sequence

for ee = 1:K
    sigma = sqrt(((1) / (2 * SNRL(ee))));  % Complex noise variance: 2 * sigma^2
    kk = 1;
    for iteration = 1:bits_total / bpcu
        bits = bit_seq(kk:kk + bpcu - 1);

        % Generate Rayleigh fading channels.
        H11 = (randn(N, Nt) + 1i * randn(N, Nt)) / sqrt(2);  % Tx-to-RIS channel
        H22 = (randn(Nr, N) + 1i * randn(Nr, N)) / sqrt(2);  % RIS-to-Rx channel

        % Map bits to an antenna combination.
        dec = sum(bits(1:end) .* [2.^(p1 - 1:-1:0)]) + 1;

        % Random antenna-combination selection.
        rowrank = randperm(size(aTx, 1));
        AS_index1 = aTx(rowrank, :);
        AS_index = AS_index1(1:Nc, :);

        % Assign RIS elements to the selected antennas.
        index = AS_index(dec, :);
        Na = nnz(index);
        index(find(index == 0)) = [];
        h = H22(index, :);
        G2 = zeros(1, N);
        for t = 1:N
            G2(1, t) = h(floor(mod(t, Na) + 1), t);
        end

        % RIS phase adjustment.
        Phi = exp(-1i * angle(G2));

        % Received signal.
        Es = 1;
        n_B = (sigma) .* (randn(Nr, Nt) + 1i * randn(Nr, Nt));
        A1 = H22 .* Phi * H11;
        y_B = sqrt(Es) * A1 + n_B;

        % ML detection.
        metrics = zeros(1, Nc);
        for count1 = 1:Nc
            index_1 = AS_index(count1, :);
            Na_1 = nnz(index_1);
            index_1(find(index_1 == 0)) = [];
            h_1 = H22(index_1, :);
            G2_1 = zeros(1, N);
            for t = 1:N
                G2_1(1, t) = h_1(floor(mod(t, Na_1) + 1), t);
            end
            Phi1 = exp(-1i * angle(G2_1));
            A2 = H22 .* Phi1 * H11;
            metrics(count1) = norm((y_B - sqrt(Es) * A2), 2);
        end

        aa1 = min(metrics);
        [a1] = ind2sub(size(metrics), find(metrics(:) == aa1));

        % Recover antenna-combination bits (MSB first).
        dec1 = zeros(1, p1);
        for pp = 1:p1
            dec1(pp) = mod(floor((a1 - 1) / (2^(p1 - pp))), 2);
        end

        decoded(kk:kk + p1 - 1) = dec1;
        kk = kk + p1;
    end

    % Estimate BER.
    bit_errors = sum(xor(bit_seq, decoded));
    BER(ee) = bit_errors / bits_total;
    disp(BER);
end

%% BER Results
toc;
semilogy(SNR, BER, '-gp', 'LineWidth', 1);
grid on;
xlabel('SNR-dB');
ylabel('BER');
hold on;
