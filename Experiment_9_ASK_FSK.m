%% EXPERIMENT 9 - ASK AND FSK
% Digital Communication
% ASK and BFSK generation, coherent/noncoherent detection,
% tone-spacing validation and BER analysis

clc;
clear;
close all;

rng(10);                         % Reproducible results

%% ==============================================================
% 1. BASIC PARAMETERS
% ==============================================================

Nbits = 1000;                    % Number of information bits
Tb = 1;                          % Bit duration
Fs = 1000;                       % Sampling frequency
Ns = Fs*Tb;                      % Samples per bit

fc = 20;                         % ASK carrier frequency
f1 = 10;                         % BFSK tone 1
f2 = 11;                         % BFSK tone 2

A = 1;                           % Signal amplitude

bits = randi([0 1],1,Nbits);

t_bit = (0:Ns-1)/Fs;

%% ==============================================================
% 2. ASK PASSBAND TRANSMITTER
% ==============================================================

ask = zeros(1,Nbits*Ns);

for k = 1:Nbits
    index = (k-1)*Ns + (1:Ns);
    
    if bits(k) == 1
        ask(index) = A*cos(2*pi*fc*t_bit);
    else
        ask(index) = 0;
    end
end

t = (0:length(ask)-1)/Fs;

%% ==============================================================
% 3. BFSK PASSBAND TRANSMITTER
% ==============================================================

bfsk = zeros(1,Nbits*Ns);

for k = 1:Nbits
    index = (k-1)*Ns + (1:Ns);
    
    if bits(k) == 0
        bfsk(index) = A*cos(2*pi*f1*t_bit);
    else
        bfsk(index) = A*cos(2*pi*f2*t_bit);
    end
end

%% ==============================================================
% 4. PASSBAND WAVEFORMS
% ==============================================================

figure;

subplot(2,1,1);
plot(t(1:10*Ns),ask(1:10*Ns),'LineWidth',1);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('ASK Passband Waveform');

subplot(2,1,2);
plot(t(1:10*Ns),bfsk(1:10*Ns),'LineWidth',1);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('BFSK Passband Waveform');

%% ==============================================================
% 5. POWER SPECTRA
% ==============================================================

L = length(ask);
f = (-L/2:L/2-1)*(Fs/L);

ASK_FFT = fftshift(fft(ask));
BFSK_FFT = fftshift(fft(bfsk));

ASK_PSD = abs(ASK_FFT).^2/L;
BFSK_PSD = abs(BFSK_FFT).^2/L;

figure;

subplot(2,1,1);
plot(f,ASK_PSD,'LineWidth',1);
grid on;
xlim([-50 50]);
xlabel('Frequency (Hz)');
ylabel('Power');
title('ASK Power Spectrum');

subplot(2,1,2);
plot(f,BFSK_PSD,'LineWidth',1);
grid on;
xlim([-50 50]);
xlabel('Frequency (Hz)');
ylabel('Power');
title('BFSK Power Spectrum');

%% ==============================================================
% 6. COHERENT ASK DETECTION
% ==============================================================

EbN0_dB = 10;

% Signal energy per bit
Es_ASK = sum((A*cos(2*pi*fc*t_bit)).^2)/Fs;

EbN0 = 10^(EbN0_dB/10);
N0 = Es_ASK/EbN0;

noise_std = sqrt(N0*Fs/2);

ask_rx = ask + noise_std*randn(size(ask));

ask_corr = zeros(1,Nbits);
ask_detected = zeros(1,Nbits);

carrier = cos(2*pi*fc*t_bit);

for k = 1:Nbits
    index = (k-1)*Ns + (1:Ns);
    
    ask_corr(k) = sum(ask_rx(index).*carrier)/Ns;
    
    if ask_corr(k) > A/4
        ask_detected(k) = 1;
    else
        ask_detected(k) = 0;
    end
end

%% ==============================================================
% 7. ASK CORRELATOR OUTPUT
% ==============================================================

figure;

stem(1:50,ask_corr(1:50),'filled');
grid on;
xlabel('Bit Number');
ylabel('Correlator Output');
title('ASK Coherent Correlator Output');

%% ==============================================================
% 8. ASK DECISION-STATISTIC HISTOGRAM
% ==============================================================

figure;

histogram(ask_corr(bits==0),30,'Normalization','pdf');
hold on;
histogram(ask_corr(bits==1),30,'Normalization','pdf');
grid on;
xlabel('Decision Statistic');
ylabel('Probability Density');
title('ASK Decision-Statistic Histogram');
legend('Bit 0','Bit 1');

%% ==============================================================
% 9. ASK BER
% ==============================================================

ASK_BER = mean(bits ~= ask_detected);

fprintf('\n====================================================\n');
fprintf('ASK RESULTS\n');
fprintf('====================================================\n');
fprintf('Eb/N0 = %.1f dB\n',EbN0_dB);
fprintf('ASK BER = %.6f\n',ASK_BER);

%% ==============================================================
% 10. BFSK BASIS SIGNALS
% ==============================================================

phi1 = cos(2*pi*f1*t_bit);
phi2 = cos(2*pi*f2*t_bit);

% Normalize the basis signals
phi1 = phi1/sqrt(sum(phi1.^2));
phi2 = phi2/sqrt(sum(phi2.^2));

inner_product = sum(phi1.*phi2);

fprintf('\n====================================================\n');
fprintf('MANDATORY BFSK INNER-PRODUCT VALIDATION\n');
fprintf('====================================================\n');
fprintf('f1 = %.2f Hz\n',f1);
fprintf('f2 = %.2f Hz\n',f2);
fprintf('Tone spacing = %.2f Hz\n',abs(f2-f1));
fprintf('Delta-f * Tb = %.4f\n',abs(f2-f1)*Tb);
fprintf('Inner product = %.8f\n',inner_product);

if abs(inner_product) < 0.01
    fprintf('Result: The BFSK basis signals are approximately ORTHOGONAL.\n');
else
    fprintf('Result: The BFSK basis signals are NOT orthogonal for this spacing.\n');
end

%% ==============================================================
% 11. COHERENT BFSK DETECTION
% ==============================================================

EbN0_dB = 10;
EbN0 = 10^(EbN0_dB/10);

Eb_BFSK = A^2*Tb/2;
N0 = Eb_BFSK/EbN0;

noise_std = sqrt(N0*Fs/2);

bfsk_rx = bfsk + noise_std*randn(size(bfsk));

corr1 = zeros(1,Nbits);
corr2 = zeros(1,Nbits);
bfsk_coherent = zeros(1,Nbits);

for k = 1:Nbits
    
    index = (k-1)*Ns + (1:Ns);
    
    corr1(k) = sum(bfsk_rx(index).*cos(2*pi*f1*t_bit))/Ns;
    corr2(k) = sum(bfsk_rx(index).*cos(2*pi*f2*t_bit))/Ns;
    
    if corr2(k) > corr1(k)
        bfsk_coherent(k) = 1;
    else
        bfsk_coherent(k) = 0;
    end
end

%% ==============================================================
% 12. BFSK CORRELATOR OUTPUTS
% ==============================================================

figure;

plot(1:50,corr1(1:50),'o-','LineWidth',1);
hold on;
plot(1:50,corr2(1:50),'s-','LineWidth',1);
grid on;
xlabel('Bit Number');
ylabel('Correlator Output');
title('BFSK Coherent Correlator Outputs');
legend('f_1 Correlator','f_2 Correlator');

%% ==============================================================
% 13. BFSK NONCOHERENT ENERGY DETECTOR
% ==============================================================

energy1 = zeros(1,Nbits);
energy2 = zeros(1,Nbits);
bfsk_noncoherent = zeros(1,Nbits);

for k = 1:Nbits
    
    index = (k-1)*Ns + (1:Ns);
    
    r1 = bfsk_rx(index).*cos(2*pi*f1*t_bit);
    r2 = bfsk_rx(index).*cos(2*pi*f2*t_bit);
    
    energy1(k) = sum(r1.^2);
    energy2(k) = sum(r2.^2);
    
    if energy2(k) > energy1(k)
        bfsk_noncoherent(k) = 1;
    else
        bfsk_noncoherent(k) = 0;
    end
end

%% ==============================================================
% 14. BFSK DECISION-STATISTIC HISTOGRAM
% ==============================================================

bfsk_stat = energy2 - energy1;

figure;

histogram(bfsk_stat(bits==0),30,'Normalization','pdf');
hold on;
histogram(bfsk_stat(bits==1),30,'Normalization','pdf');
grid on;
xlabel('Decision Statistic: E_2 - E_1');
ylabel('Probability Density');
title('BFSK Noncoherent Decision-Statistic Histogram');
legend('Bit 0','Bit 1');

%% ==============================================================
% 15. BFSK BER
% ==============================================================

BER_BFSK_COH = mean(bits ~= bfsk_coherent);
BER_BFSK_NONCOH = mean(bits ~= bfsk_noncoherent);

fprintf('\n====================================================\n');
fprintf('BFSK RESULTS\n');
fprintf('====================================================\n');
fprintf('Eb/N0 = %.1f dB\n',EbN0_dB);
fprintf('Coherent BFSK BER = %.6f\n',BER_BFSK_COH);
fprintf('Noncoherent BFSK BER = %.6f\n',BER_BFSK_NONCOH);

%% ==============================================================
% 16. BER CURVES
% ==============================================================

EbN0_range = 0:2:14;

BER_ASK = zeros(size(EbN0_range));
BER_BFSK_COH = zeros(size(EbN0_range));
BER_BFSK_NONCOH = zeros(size(EbN0_range));

Nber = 20000;

bits_ber = randi([0 1],1,Nber);

ask_clean = zeros(1,Nber*Ns);
bfsk_clean = zeros(1,Nber*Ns);

for k = 1:Nber
    
    index = (k-1)*Ns + (1:Ns);
    
    if bits_ber(k) == 1
        ask_clean(index) = A*cos(2*pi*fc*t_bit);
        bfsk_clean(index) = A*cos(2*pi*f2*t_bit);
    else
        ask_clean(index) = 0;
        bfsk_clean(index) = A*cos(2*pi*f1*t_bit);
    end
end

for m = 1:length(EbN0_range)
    
    EbN0 = 10^(EbN0_range(m)/10);
    
    %% ASK
    
    Eb_ASK = A^2*Tb/4;
    N0 = Eb_ASK/EbN0;
    sigma = sqrt(N0*Fs/2);
    
    rx = ask_clean + sigma*randn(size(ask_clean));
    
    stat = zeros(1,Nber);
    
    for k = 1:Nber
        
        index = (k-1)*Ns + (1:Ns);
        stat(k) = sum(rx(index).*carrier)/Ns;
        
    end
    
    detected = stat > A/4;
    BER_ASK(m) = mean(bits_ber ~= detected);
    
    %% BFSK
    
    Eb_BFSK = A^2*Tb/2;
    N0 = Eb_BFSK/EbN0;
    sigma = sqrt(N0*Fs/2);
    
    rx = bfsk_clean + sigma*randn(size(bfsk_clean));
    
    stat1 = zeros(1,Nber);
    stat2 = zeros(1,Nber);
    
    for k = 1:Nber
        
        index = (k-1)*Ns + (1:Ns);
        
        stat1(k) = sum(rx(index).*cos(2*pi*f1*t_bit))/Ns;
        stat2(k) = sum(rx(index).*cos(2*pi*f2*t_bit))/Ns;
        
    end
    
    detected = stat2 > stat1;
    BER_BFSK_COH(m) = mean(bits_ber ~= detected);
    
    %% Noncoherent BFSK
    
    e1 = zeros(1,Nber);
    e2 = zeros(1,Nber);
    
    for k = 1:Nber
        
        index = (k-1)*Ns + (1:Ns);
        
        z1 = rx(index).*cos(2*pi*f1*t_bit);
        z2 = rx(index).*cos(2*pi*f2*t_bit);
        
        e1(k) = sum(z1.^2);
        e2(k) = sum(z2.^2);
        
    end
    
    detected = e2 > e1;
    BER_BFSK_NONCOH(m) = mean(bits_ber ~= detected);
    
end

%% ==============================================================
% 17. BER CURVE
% ==============================================================

figure;

semilogy(EbN0_range,BER_ASK,'o-','LineWidth',1.5);
hold on;
semilogy(EbN0_range,BER_BFSK_COH,'s-','LineWidth',1.5);
semilogy(EbN0_range,BER_BFSK_NONCOH,'^-','LineWidth',1.5);

grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance of ASK and BFSK');
legend('ASK Coherent','BFSK Coherent','BFSK Noncoherent');
ylim([1e-5 1]);

%% ==============================================================
% 18. TONE-SPACING VARIATION
% ==============================================================

spacing_values = [0.5 1 1.5 2 3];
spacing_inner_product = zeros(size(spacing_values));

fprintf('\n====================================================\n');
fprintf('TONE-SPACING VALIDATION\n');
fprintf('====================================================\n');

fprintf('Before variation: Increasing tone spacing should reduce\n');
fprintf('correlation between the BFSK basis signals. Orthogonality\n');
fprintf('is expected when Delta-f * Tb is an integer.\n\n');

for m = 1:length(spacing_values)
    
    test_f2 = f1 + spacing_values(m);
    
    test_phi2 = cos(2*pi*test_f2*t_bit);
    test_phi2 = test_phi2/sqrt(sum(test_phi2.^2));
    
    spacing_inner_product(m) = sum(phi1.*test_phi2);
    
    fprintf('Spacing = %.2f Hz, Delta-f*Tb = %.2f, Inner product = %.8f\n',...
        spacing_values(m),spacing_values(m)*Tb,...
        spacing_inner_product(m));
end

figure;

plot(spacing_values,abs(spacing_inner_product),'o-','LineWidth',1.5);
grid on;
xlabel('Tone Spacing, Delta f (Hz)');
ylabel('|Inner Product|');
title('BFSK Orthogonality Validation');

%% ==============================================================
% 19. FINAL OBSERVATION
% ==============================================================

fprintf('\n====================================================\n');
fprintf('OBSERVATION AND INTERPRETATION\n');
fprintf('====================================================\n');

fprintf('1. ASK: Increasing Eb/N0 should reduce BER because noise\n');
fprintf('   becomes less significant relative to the signal.\n');

fprintf('2. BFSK: Coherent detection uses phase-synchronized\n');
fprintf('   correlators, while noncoherent detection uses energy.\n');

fprintf('3. BFSK orthogonality requires Delta-f * Tb to be an\n');
fprintf('   integer for rectangular tone pulses.\n');

fprintf('4. The inner-product test above is used to diagnose whether\n');
fprintf('   the selected BFSK tones are orthogonal.\n');

fprintf('5. Any non-zero inner product indicates imperfect\n');
fprintf('   orthogonality and can affect coherent detection.\n');

fprintf('\nExperiment 9 simulation completed successfully.\n');