%%---File Setup And Main Function
clc;
clear;
close all;

% --- Convert to WAV if needed ---
input_file = 'Amr Diab - Lama Abeltak_ عمرو دياب - لما قابلتك(M4A_128K).wav';
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
start_sample = fs_original * 1+ 1;        
end_sample = min(length(Signal), fs_original * 20);  
Signal = Signal(start_sample:end_sample);

% --- User Inputs ---
fs = input('Enter the desired sampling frequency (e.g., 40000): ');
L = input('Enter the number of quantization levels (e.g., 4, 8, 64): ');

% --- Sampling ---
t = (0:length(Signal)-1)/fs_original;
[t_s, Signal_s] = sampler(Signal, t, fs);


% --- Plot Sampler Output (Zoomed to First 20s) ---
plot_duration = 20;
idx_full = t <= plot_duration;
idx_sampled = t_s <= plot_duration;

figure;
subplot(2, 1, 1);
plot(t(idx_full), Signal(idx_full), 'k');
xlabel('Time (s)'); ylabel('Amplitude'); title('Original Signal'); grid on;

subplot(2, 1, 2);
stem(t_s(idx_sampled), Signal_s(idx_sampled), 'b', 'filled');
xlabel('Time (s)'); ylabel('Amplitude'); title('Sampled Signal'); grid on;

% --- Quantization ---
[Quantized_input, quant_levels, mse, bitstream, mp_min, mp_max] = quantizer(Signal_s, t_s, L, []);

% --- Encoding ---
[waveform, t_pulse, Tb, amplitude, type] = encoder(bitstream, fs, L);

% --- Show first 10 bits waveform ---
figure;
stairs(t_pulse(1:10), waveform(1:10), 'LineWidth', 2);
title('First 10 Bits of NRZ-Encoded Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

% --- Decoding (using original decoder) ---
x_dec = decoder(waveform, quant_levels, amplitude, type);

% --- Plot quantized vs decoded ---
figure;
plot(t_s, Quantized_input, 'b', t_s, x_dec, 'r--');
legend('Quantized Signal', 'Decoded Signal');
title('Quantized vs Decoded Output'); grid on;

% --- Save output Signal ---
output_name = sprintf('decoded_user_fs%d_L%d.wav', fs, L);
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


function [waveform, t_pulse, Tb, amplitude, type] = encoder(bitstream, fs, L)
    bits_per_sample = ceil(log2(L));
    default_Tb = 1 / (fs * bits_per_sample);

    fprintf('Auto-calculated bit duration Tb = %.6f s\n', default_Tb);
    use_default = input('Use this value? (y/n): ', 's');

    if strcmpi(use_default, 'y')
        Tb = default_Tb;
    else
        Tb = input('Enter your custom NRZ bit duration (in seconds): ');
    end

    amplitude = input('Enter NRZ pulse amplitude: ');
    type = input('Choose NRZ type (unipolar / polar): ', 's');

    t_pulse = 0:Tb:(length(bitstream) * Tb - Tb);
    waveform = zeros(size(t_pulse));

    for i = 1:length(bitstream)
        if strcmpi(type, 'unipolar')
            waveform(i) = amplitude * bitstream(i);
        elseif strcmpi(type, 'polar')
            waveform(i) = amplitude * (2 * bitstream(i) - 1);
        else
            error('Invalid NRZ type. Choose "unipolar" or "polar".');
        end
    end
end
%%
%--Decoder Function Definition--%
function x_decoded = decoder(waveform, quant_levels, amplitude, type)
    bits_per_sample = ceil(log2(length(quant_levels)));
    bitstream = zeros(1, length(waveform));

    for i = 1:length(waveform)
        if strcmpi(type, 'unipolar')
            bitstream(i) = round(waveform(i) / amplitude);
        elseif strcmpi(type, 'polar')
            bitstream(i) = waveform(i) > 0;
        else
            error('Invalid NRZ type.');
        end
    end

    n_samples = length(bitstream) / bits_per_sample;
    x_decoded = zeros(1, n_samples);
    for i = 1:n_samples
        bits = bitstream((i-1)*bits_per_sample + 1 : i*bits_per_sample);
        idx = bin2dec(char(bits + '0')) + 1;
        x_decoded(i) = quant_levels(idx);
    end

    t = linspace(0, length(x_decoded)/n_samples, length(x_decoded));
    figure;
    stem(t, x_decoded, 'm');
    title('Decoded Quantized Signal');
    xlabel('Time [s]'); ylabel('Amplitude'); grid on;
end
