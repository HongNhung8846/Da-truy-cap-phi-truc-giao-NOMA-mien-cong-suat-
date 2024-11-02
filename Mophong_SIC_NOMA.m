clc; clear variables; close all;
N = 10^5;               %Number of monte carlo simulations
SNR = 0:30;             %SNR range in dB
snr = db2pow(SNR);      %SNR range in linear scale
%Generate random data bits for transmission
x1 = randi([0 1],1,N);  %Data bits of user 1
x2 = randi([0 1],1,N);  %Data bits of user 2
x3 = randi([0 1],1,N);  %Data bits of user 3
%Do BPSK modulation of data
xmod1 = 2*x1 - 1;
xmod2 = 2*x2 - 1;
xmod3 = 2*x3 - 1;
%Set power weights for users
a1 = 0.6; a2 = 0.25; a3 = 0.15;
%Do superposition coding
x = sqrt(a1)*xmod1 + sqrt(a2)*xmod2 + sqrt(a3)*xmod3
%Add AWGN to x (Transmit x through an AWGN channel)
for u = 1:length(snr)
    y1 = awgn(x,SNR(u),'measured');  %Received signal at user 1 corrupted by AWGN
    y2 = awgn(x,SNR(u),'measured'); %Received signal at user 2 corrupted by AWGN
    y3 = awgn(x,SNR(u),'measured'); %Received signal at user 2 corrupted by AWGN  
    %AT USER 1
    %Direct decoding of x from y1
    x1_hat = ones(1,N); %just a buffer 
    x1_hat(y1 < 0) = 0;  
    %AT USER 2
    %Direct decoding of x from y2
    x11_est = ones(1,N); %just a buffer 
    x11_est(y2 < 0) = 0; %Estimate user 1's signal first
    x11_est(x11_est == 0) = -1; %Remodulate user 1's signal
    %Subtract remodulated x11_est component from y2
    rem = y2 - sqrt(a1)*x11_est;  
    %Decode x2 from rem
    x2_hat = zeros(1,N);
    x2_hat(rem>0) = 1; 
    % AT USER 3
    %Direct decoding of x from y3
    x111_est = ones(1, N);
    x111_est(y3 < 0) = 0;
    x111_est(x111_est == 0) = -1;   
    %Subtract remodulated x111_est component from y3
    rem3 = y3 - sqrt(a1) * x111_est;
    x22_est = ones(1, N);
    x22_est(rem3 < 0) = 0;
    x22_est(x22_est == 0) = -1;
    rem3_final = rem3 - sqrt(a2) * x22_est;
    %Decode x2 from rem
    x3_hat = ones(1, N);
    x3_hat(rem3_final < 0) = 0;  
    %Estimate BER
    ber1(u) = biterr(x1,x1_hat)/N;
    fprintf('ber1: %.4f\n', ber1);
    ber2(u) = biterr(x2,x2_hat)/N;
    fprintf('ber2: %.4f\n',ber2);
    ber3(u) = biterr(x3,x3_hat)/N;
    fprintf('ber3: %.4f\n',ber3);
end
%plot BER curves
figure(1);
semilogy(SNR, ber1, 'linewidth', 1.5); hold on;
semilogy(SNR, ber2, 'linewidth', 1.5); grid on;
semilogy(SNR, ber3, 'linewidth', 1.5); grid on;
legend('User 1 \alpha_1 = 0.6','User 2 \alpha_2 = 0.25', 'User 3 \alpha_3 = 0.15');
xlabel('SNR (dB)');
ylabel('BER');
title('BER graph for NOMA in AWGN channel');
 
figure(2);
subplot(3,2,1);
plot(x(1:1000), 'LineWidth', 1.5); title('Original Transmitted Signal'); xlabel('Sample'); ylabel('Amplitude');
subplot(3, 2, 2);
plot(x1_hat(1:1000), 'LineWidth', 1.5); title('Decoded Signal at User 1'); xlabel('Sample'); ylabel('Bit Value');
subplot(3,2,3);
plot(rem(1:1000), 'LineWidth', 1.5); title('Remaining Signal after Removing User 1 at User 2'); xlabel('Sample'); ylabel('Amplitude');
subplot(3,2,4);
plot(x2_hat(1:1000), 'LineWidth', 1.5); title('Decoded Signal at User 2'); xlabel('Sample'); ylabel('Bit Value');
subplot(3,2,5);
plot(rem3_final(1:1000), 'LineWidth', 1.5); title('Remaining Signal after Removing User 1 and User 2 at User 3'); xlabel('Sample'); ylabel('Amplitude');
subplot(3,2,6);
plot(x3_hat(1:1000), 'LineWidth', 1.5); title('Decoded Signal at User 3'); xlabel('Sample'); ylabel('Bit Value');
