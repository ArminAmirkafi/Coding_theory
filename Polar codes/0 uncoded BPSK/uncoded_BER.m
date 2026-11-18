clc; clear; close all;

% Simulation Parameters
N = 256;                  
K = 128;                 
EbNodB_range = 0:0.5:6;   
max_frames = 20000;       
max_errors = 100;         

disp('Running Uncoded BPSK...');
BER_uncoded = zeros(1, length(EbNodB_range));

for idx = 1:length(EbNodB_range)
    EbNodB = EbNodB_range(idx);

    snr_linear = 10^(EbNodB / 10);
    sigma = sqrt(1 / (2 * snr_linear));
    
    bit_errors = 0;
    frames = 0;
    
    while bit_errors < max_errors && frames < max_frames
        msg = randi([0 1], 1, K);
        
        % BPSK 
        s = 1 - 2 * msg;
        
        %  AWGN
        noise = sigma * randn(1, K);
        r = s + noise; 
        
        % 4. Hard Demodulation: r < 0 maps back to 1, r >= 0 maps to 0
        msg_hat = r < 0;
        
        % Count errors
        bit_errors = bit_errors + sum(msg ~= msg_hat);
        frames = frames + 1;
    end
    
    BER_uncoded(idx) = bit_errors / (frames * K);
    fprintf('Uncoded | Eb/No = %4.1f dB | BER = %9.2e | Frames = %d\n', EbNodB, BER_uncoded(idx), frames);
end

% Plotting
figure;
% Plot Uncoded
semilogy(EbNodB_range, BER_uncoded, 'k-s', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'DisplayName', 'Uncoded BPSK');
hold on;

% Also plot theoretical Uncoded BPSK for sanity check
theoretical_ber = 0.5 * erfc(sqrt(10.^(EbNodB_range/10)));
semilogy(EbNodB_range, theoretical_ber, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Theoretical BPSK');

grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance: Uncoded vs Theoretical BPSK');
legend('show', 'Location', 'southwest');
ylim([1e-5 1]); % Lock y-axis to standard BER viewing range