



clear

%% Test signal

fs = 16000;
T = 1/fs;
t1_vct = (0 : 1/fs : 10 - 1/fs)';

x1 = 3 * sin (2 * pi * 10 * t1_vct);
x2 = 2 * sin (2 * pi * 24 * t1_vct);
x3 = 1 * randn([numel(t1_vct), 1]);

x = [x1; x2; x3];
n = numel(x);


% Power using its PSD from rFFT
psd_rfft_b = rfft_psd(x, fs, [], 'blackmanharris');
f_step = psd_rfft_b.freq_axis(2);
power_psd_rfft_x_bh = f_step * sum(psd_rfft_b.PSD);
figure()
plot_psd_data(psd_rfft_b)

%% Power from STFFT Spectrogram (Hamming window)
w_size =  1 * fs;
w_shift = 0.5 * w_size;
win_funct = 'hamming';
rfft_spect_h = strfft_spectrogram(x, fs, w_size, w_shift, [], win_funct );
power_spect_h = sum(sum(rfft_spect_h.power_spectrogram)) * rfft_spect_h.freq_delta * rfft_spect_h.time_delta;
figure()
plot_spectrogram_data(rfft_spect_h);

%% Power from Modulation Spectrogram STFFT
w_size =  1 * fs;
w_shift = 0.5 * w_size;
rfft_mod_b = strfft_modulation_spectrogram(x, fs, w_size, w_shift, [], win_funct, [], win_funct);
power_mod = sum(sum(rfft_mod_b.power_modulation_spectrogram) * rfft_mod_b.freq_delta * rfft_mod_b.freq_mod_delta);
figure()
plot_modulation_spectrogram_data(rfft_mod_b)
