%% TODO

% Include filter response in relation to arrows

% Experimental results: draw graphs corresponding to transformation /
% synthesis of the data (e.g. factor loadings -> freq. boundaries ->
% curves comparing music / speech...

% Draw half-wave rectification + LP-filtering

% Draw envelope-modulated sinewave

% Compute waveform + spectrogram of the 9 digits of french and display them

% Build coloured confusion matrices


subplot(1, 1, 1); % Wipe drawing surface
set(gcf, 'Color', [1, 1, 1]);

%% Requirements
% voicebox toolbox () for filterbank decomposition
% vocoder toolbox (Etienne Gaudrain) for envelope extraction (or my own)?

%% Load audio file
file = "~/work/sounds/databases/lscp/french/frl1151.aiff";
[sig, fs] = audioread(file, 'double');

% Generate random noise for testing
%sig = randn(1000, 1);
%fs = 44100;

% Generate temporal vector
t = linspace(0, fs/length(sig), length(sig));

%%
gh = gcf();
pos = get(gh, 'Position'); %// gives x left, y bottom, width, height
fw = pos(3);
fh = pos(4);


ph = 8;
pw = 15;

%% Start plotting commands: Initial signal
%subplot(ph, pw, (i+0):(i+2));

i = 1+(3*pw); % Row number computation
subplot(ph, pw, (i+0):(i+2));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')


%% Row 1
i = 1; % Row number computation
% Signal
subplot(ph, pw, (i+4):(i+6));
plot(t, sig, 'Color', 'blue');
ylim([min(sig), max(sig)]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')

% Envelope
subplot(ph, pw, (i+8):(i+10));
srms = matlab.tall.movingWindow('rms', round(fs/100), sig);
%srms = 20*log10(srms/max(srms));
srms = 20*log10(srms/min(srms));
plot(t, srms);
%axis off;
ylim(1.*[0, max(srms)]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')

% Alternate
subplot(ph, pw, (i+12):(i+14));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')


%% Row 2
i = 1+(2*pw); % Row number computation
subplot(ph, pw, (i+4):(i+6));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')
subplot(ph, pw, (i+8):(i+10));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')
subplot(ph, pw, (i+12):(i+14));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')

%% Row 3
i = 1+(4*pw); % Row number computation
subplot(ph, pw, (i+4):(i+6));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')
subplot(ph, pw, (i+8):(i+10));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')
subplot(ph, pw, (i+12):(i+14));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')

%% Row 4
i = 1+(6*pw); % Row number computation
subplot(ph, pw, (i+4):(i+6));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')
subplot(ph, pw, (i+8):(i+10));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')
subplot(ph, pw, (i+12):(i+14));
ylim([-1, 1]);
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', 'box', 'off')

%subplot(ph, pw, 14:15);
%subplot(ph, pw, 17:18);
%subplot(ph, pw, 20:22);

%% Arrows
% Filter etc
etc = annotation('textbox', [0.325, 0.1, 0.5, 0.5], 'String', '[...]', 'FitBoxToText','on', 'FontSize', 24);
etc.EdgeColor = [1 1 1];

% Filterbank
etc = annotation('textbox', [0.375, 0.495, 0.5, 0.5], 'String', 'Filterbank', 'FitBoxToText','on', 'FontSize', 16);
etc.EdgeColor = [1 1 1];

% Envelopes
etc = annotation('textbox', [0.575, 0.495, 0.5, 0.5], 'String', 'Envelope', 'FitBoxToText','on', 'FontSize', 16);
etc.EdgeColor = [1 1 1];

% Energy information or Envelope modulation?
etc = annotation('textbox', [0.775, 0.495, 0.5, 0.5], 'String', 'Energy', 'FitBoxToText','on', 'FontSize', 16);
etc.EdgeColor = [1 1 1];

% Coordinates : [x1, x2], [y1, y2]
xpos = [0.495, 0.495+0.03];
yref = 1/ph;
yref
arh=annotation('arrow', xpos, (ph-0.9)*repmat(yref, 2, 1));
arh=annotation('arrow', xpos, (ph-2.58)*repmat(yref, 2, 1));
arh=annotation('arrow', xpos, (ph-4.27)*repmat(yref, 2, 1));
arh=annotation('arrow', xpos, (ph-5.96)*repmat(yref, 2, 1));

% Coordinates : [x1, x2], [y1, y2]
xpos = [0.705, 0.705+0.03];
yref = 1/ph;
yref
arh=annotation('arrow', xpos, (ph-0.9)*repmat(yref, 2, 1));
arh=annotation('arrow', xpos, (ph-2.58)*repmat(yref, 2, 1));
arh=annotation('arrow', xpos, (ph-4.27)*repmat(yref, 2, 1));
arh=annotation('arrow', xpos, (ph-5.96)*repmat(yref, 2, 1));

%arh.LineWidth = 1;
%arh.Color = [1,0,1];

%arh=annotation('arrow', [0.3,0.35],[0.85,0.85] )
%arh.LineWidth = 1;
%arh.Color = [1,0,1];
