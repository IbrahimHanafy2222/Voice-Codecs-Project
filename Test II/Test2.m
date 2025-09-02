%%---File Setup And Main Function
clc;
clear;
close all;

% --- Convert to WAV if needed ---
input_file = 'Amr Diab - Lama Abeltak_ عمرو دياب - لما قابلتك(M4A_128K).m4a';
[~, name, ~] = fileparts(input_file);
converted_file = [name, '.wav'];

if ~isfile(converted_file)
    [y, fs_conv] = audioread(input_file);
    audiowrite(converted_file, y, fs_conv);
    fprintf('Converted "%s" to "%s"\n', input_file, converted_file);
end

% --- Load and preprocess audio ---
[Signal, fs_original] = audioread(converted_file);
Signal = mean(Signal, 2);  % Convert stereo to mono
start_sample = fs_original * 1 + 1;
end_sample = min(length(Signal), fs_original * 20);
Signal = Signal(start_sample:end_sample);
sound(Signal,44100)

% --- Parameters ---
fs = 40000;
L = 64;
                % samples per bit
pulse_amplitude = 2;
N0 = 1;                   % Noise variance

% --- Sampling ---
t = (0:length(Signal)-1)/fs_original;
[t_s, Signal_s] = sampler(Signal, t, fs);

% --- Plot Sampled Signal ---
figure;
subplot(2, 1, 1);
plot(t(t <= 20), Signal(t <= 20), 'k');
xlabel('Time (s)'); ylabel('Amplitude'); title('Original Signal'); grid on;

subplot(2, 1, 2);
stem(t_s(t_s <= 20), Signal_s(t_s <= 20), 'b', 'filled');
xlabel('Time (s)'); ylabel('Amplitude'); title('Sampled Signal'); grid on;

% --- Quantization ---
[Quantized_input, quant_levels, mse, bitstream, mp_min, mp_max] = quantizer(Signal_s, t_s, L, []);

% --- Encoding ---
[waveform, t_pulse, Tb, amplitude, type,n] = encoder(bitstream, fs, L);

% --- Channel (AWGN) ---
noisy_signal = AWGN_channel(t_pulse, waveform, N0);

% --- Regenerative Repeater ---
regenerated_PCM_signal = regenerative_repeater_NRZ(t_pulse, noisy_signal, n, pulse_amplitude);

% --- Decoder ---
x_dec = decoder(regenerated_PCM_signal, quant_levels, amplitude, type, n);

% --- Plot decoded vs quantized ---
t_dec = linspace(0, t_s(end), length(x_dec));
figure;
plot(t_s, Quantized_input, 'b', t_dec, x_dec, 'r--');
legend('Quantized Signal', 'Decoded Signal');
title('Quantized vs Decoded Output');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% --- Save output ---
output_name = sprintf('unipolar_1_decoded_user_fs%d_L%d.wav', fs, L);
audiowrite(output_name, x_dec, fs);
fprintf('Decoded audio saved as "%s"\n', output_name);

%%
%--Sampler--%
function [t_sampled, Signal_sampled] = sampler(signal, t_original, fs_target)
    duration = t_original(end);
    t_sampled = 0:1/fs_target:duration;
    Signal_sampled = interp1(t_original, signal, t_sampled, 'linear');
end
%%
%--Quantizer Function Definition--% 
function [Quantized_input, quant_levels, mse, bitstream, mp_min, mp_max] = quantizer(Signal, t, L, ~)
    mode = input('Choose quantizer type (mid-rise / mid-tread): ', 's');
    mp_min = min(Signal);
    mp_max = max(Signal);
    delta = (mp_max - mp_min) / L;

    switch lower(mode)
        case 'mid-rise'
            min_level = mp_min + delta/2;
            max_level = mp_max - delta/2;
            quant_levels = min_level : delta : max_level;
        case 'mid-tread'
            min_level1 = mp_min + delta;
            max_level1 = mp_max;
            quant_levels = min_level1 : delta : max_level1;
        otherwise
            error('Invalid quantizer type. Choose either "mid-rise" or "mid-tread".');
    end

    Quantized_input = zeros(size(Signal));
    indices = zeros(size(Signal));
    for i = 1:length(Signal)
        [~, idx] = min(abs(Signal(i) - quant_levels));
        Quantized_input(i) = quant_levels(idx);
        indices(i) = idx - 1;
    end

    mse = mean((Signal - Quantized_input).^2);

    bits_per_sample = ceil(log2(length(quant_levels)));
    bin_matrix = dec2bin(indices, bits_per_sample) - '0';
    bitstream = reshape(bin_matrix.', 1, []);

    figure;
    plot(t, Signal, 'b', 'DisplayName', 'Input Signal'); hold on;
    stairs(t, Quantized_input, 'r', 'DisplayName', 'Quantized Signal');
    xlabel('t [sec]'); ylabel('Amplitude');
    title('Input Signal vs. Quantized Signal'); legend('show'); grid on;
end
%%
%--Encoder Function Definition--%
function [waveform, t_pulse, Tb, amplitude, type, n] = encoder(bitstream, fs, L)
    bits_per_sample = ceil(log2(L));
    default_Tb = 1 / (fs * bits_per_sample);

    fprintf('Auto-calculated bit duration Tb = %.6f s\n', default_Tb);
    use_default = input('Use this value? (y/n): ', 's');

    if strcmpi(use_default, 'y')
        Tb = default_Tb;
    else
        Tb = input('Enter your custom NRZ bit duration (in seconds): ');
    end

    amplitude = 2;
    type = input('Choose NRZ type (unipolar / polar): ', 's');

    n = 100; % number of samples per bit (oversampling factor)
    waveform = zeros(1, n * length(bitstream));

    for i = 1:length(bitstream)
        if strcmpi(type, 'unipolar')
            value = amplitude * bitstream(i);
        elseif strcmpi(type, 'polar')
            value = amplitude * (2 * bitstream(i) - 1);
        else
            error('Invalid NRZ type. Choose "unipolar" or "polar".');
        end
        waveform((i-1)*n + 1 : i*n) = value;
    end

    t_pulse = 0:Tb/n:Tb*length(bitstream)-Tb/n;

    % === Plot the Bitstream Waveform ===
    figure;
    plot(t_pulse(1:20*n), waveform(1:20*n), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['NRZ Encoded Bitstream - ', upper(type)]);
    grid on;
end

%%
%--Decoder Function Definition--%
function x_decoded = decoder(waveform, quant_levels, amplitude, type, n)
    % Number of bits used per quantized sample
    bits_per_sample = ceil(log2(length(quant_levels)));

    % Step 1: Extract bitstream from waveform using averaging
    samples_per_bit = n;
    num_bits = length(waveform) / samples_per_bit;
    bitstream = zeros(1, num_bits);

    for i = 1:num_bits
        segment = waveform((i-1)*samples_per_bit + 1 : i*samples_per_bit);
        if strcmpi(type, 'polar')
            bitstream(i) = mean(segment) > 0;
        elseif strcmpi(type, 'unipolar')
            bitstream(i) = mean(segment) > amplitude / 2;
        else
            error('Invalid NRZ type. Choose "polar" or "unipolar".');
        end
    end

    % Step 2: Convert bitstream to quantized indices
    n_samples = length(bitstream) / bits_per_sample;
    x_decoded = zeros(1, n_samples);

    for i = 1:n_samples
        bits = bitstream((i-1)*bits_per_sample + 1 : i*bits_per_sample);
        idx = bin2dec(char(bits + '0')) + 1;
        x_decoded(i) = quant_levels(idx);
    end

    % Step 3: Plot decoded signal
    t = linspace(0, length(x_decoded)/n_samples, length(x_decoded));
    figure;
    plot(t, x_decoded, 'm');
    title('Decoded Quantized Signal');
    xlabel('Time [s]');
    ylabel('Amplitude');
    grid on;
end
%% 
%--Noise Generator--%
function [noisy_signal] = AWGN_channel(t, signal, N0)
    signal_power = mean(signal.^2);
    snr_linear = signal_power / N0;

    noisy_signal = signal + sqrt(N0) * randn(size(signal));

    figure;
    plot(t(1:100*20), noisy_signal(1:100*20));   % Show 20 bits worth if n=100
    xlabel('t [sec]');
    ylabel('Amplitude');
    title('Noisy PCM signal first 20 bits');
    legend('Channel output');
    grid on;
end

%%
%--Repeater%--%
function [regenerated_PCM_signal] = regenerative_repeater_NRZ(t, PCM_signal, n, pulse_amplitude)

    A = pulse_amplitude;
    regenerated_PCM_signal = zeros(1, length(PCM_signal));
    prompt = 'Enter NRZ line code type (polar or unipolar): ';
    line_code = lower(input(prompt, 's'));

    if strcmp(line_code, 'polar')
        for i = 1:n:length(PCM_signal)-n+1
            segment = PCM_signal(i:i+n-1);
            corr_pos = dot(segment, A * ones(1, n));
            corr_neg = dot(segment, -A * ones(1, n));
            
            if corr_pos > corr_neg
                regenerated_PCM_signal(i:i+n-1) = A;
            else
                regenerated_PCM_signal(i:i+n-1) = -A;
            end
        end

    % elseif strcmp(line_code, 'unipolar')
    %     for i = 1:n:length(PCM_signal)-n+1
    %         segment = PCM_signal(i:i+n-1);
    %         corr_1 = dot(segment, A * ones(1, n));
    %         corr_0 = dot(segment, zeros(1, n));
    % 
    %         if corr_1 > corr_0 + 0.1  % small margin to account for noise
    %             regenerated_PCM_signal(i:i+n-1) = A;
    %         else
    %             regenerated_PCM_signal(i:i+n-1) = 0;
    %         end
    elseif strcmp(line_code, 'unipolar')
    for i = 1:n:length(PCM_signal)-n+1
        segment = PCM_signal(i:i+n-1);
        average = mean(segment);
        
        if average > A / 2   % Threshold set to half amplitude
            regenerated_PCM_signal(i:i+n-1) = A;
        else
            regenerated_PCM_signal(i:i+n-1) = 0;
        end

        end

    else
        error('Invalid line code type. Use "polar" or "unipolar".');
    end

    % Plotting
    nexttile
    plot(t(1:min(end, 20*n)), regenerated_PCM_signal(1:min(end, 20*n)));
    xlabel('t [sec]');
    ylabel('Amplitude');
    title(['Regenerated PCM Signal (First 20 Bits) - ', upper(line_code), ' NRZ']);
    legend('Regenerated Signal');
    grid on;
end

