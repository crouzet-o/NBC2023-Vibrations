

gaudrain_path = '~/work/src/github/others/vocoder/';
addpath(gaudrain_path);

vocbands = false;
singleorig = false;
recombinedorig = false;
carriers = false;
envs = false;

% Read a sound file



% Text variables have to be stored in cells in order to be functional
ctypes = {'noise', 'sin', 'low_noise', 'phsc'};
%ctype = char(ctypes(1)); % noise
%ctype = char(ctypes(2)); % sin
ctype = char(ctypes(4)); % phsc

HFlimit = 16000;
% Use up-scale translation (2 ERBs)
parBands = [1, 2, 3, 4, 6, 8, 12, 16, 22, 32];
%parBands = [12];

%nbbands = parBands(4); % 4 channels
%nbbands = parBands(6); % 8 channels
nbbands = parBands(7); % 12 channels
nbbands = parBands(8); % 16 channels
%nbbands = parBands(9); % 22 channels
%nbbands = parBands(10); % 32 channels

%% Create a list of filenames (a cell array);
%% Choose a fixed number of bands / type of carrier (sinewave vs. noise);
%% 

%songs = {
%    "01-maurice-ravel-bolero.mp3",
%    "02-claude-debussy-prelude-a-lapres-midi-dun-faune.mp3",
%    "03-le-sacre-du-printemps-danse-des-adolescentes-vladimir-kulenovic-belgrade-philharmonic.mp3",
%    "04-dean-martin-dream-a-little-dream-of-me.mp3",
%    "05-sarah-vaughan-summertime.mp3",
%    "06-herbie-hancock-fat-albert-rotunda.mp3",
%    "07-stan-getz-joao-gilberto-desafinado.mp3",
%    "10 - Cypress Hill - Insane In The Brain.mp3",
%    "11-michael-jackson-billie-jean-official-music-video.mp3",
%    "12-michael-jackson-beat-it-official-video.mp3",
%    "13-prince-purple-rain-official-video.mp3",
%    "14-prince-kiss-official-music-video.mp3",
%    "15-queen-bohemian-rhapsody-live-at-wembley-stadium-saturday-12-july-1986.mp3",
%    "16-nirvana-smells-like-teen-spirit.mp3",
%    "17-iggy-pop-the-passenger.mp3",
%    "18-david-bowie-the-man-who-sold-the-world.mp3",
%    "10 - Stan Getz, João Gilberto - The Girl from Ipanema - extract.flac"
%}
songPath = './music/';
songs = dir([songPath, 'origins/new/*.mp3']);
numsongs = size(songs);

%	"01 - Ravel - Boléro - extract.flac",
%	"02 - Debussy - Prélude à l'après-midi d'un faune - extract.flac",
%	"03 - Stravinsky - Le Sacre du printemps: Danse des adolescentes - extract.flac",
%        "04 - Alexinho - 20-beatbox-solo-C - final - extract.flac",
%	"04 - Dream a Little Dream of Me - Dean Martin - extract.flac",
%	"05 - Summertime - Sarah Vaughan - extract.flac",
%	"06 - Fat Albert Rotunda - Herbie Hancock - extract.flac",
%	"07 - Desafinado - Stan Getz, João Gilberto - extract.flac",
%	"10 - Cypress Hill - Insane in the Brain - extract.flac",
%	"11 - Michael Jackson - Billie Jean - extract.flac",
%	"12 - Michael Jackson - Beat It - extract.flac",
%	"13 - Prince - Purple Rain - extract.flac",
%	"14 - Prince - Kiss - extract.flac",
%	"15 - Queen - Bohemian Rhapsody (Live at Wembley) - extract.flac",
%	"16 - Nirvana - Smells like Teen Spirit - extract.flac",
%	   "17 - Iggy Pop - The Passenger - extract.flac",
%	"18 - David Bowie - The Man who Sold the World - extract.flac",
%	"19 - Stan Getz, João Gilberto - The Girl from Ipanema - extract.flac"

%songs = {"01 - Ravel - Boléro - extract.flac"};

%songs = {"19 - Stan Getz, João Gilberto - The Girl from Ipanema - extract.flac"};
%inits = 0.0001; ends = 13;

% songs = {"04 - Alexinho - 20-beatbox-solo-C - final - extract.flac"};

%% 	   "08 - Darbuka and Riq Jam Istanbul 2014 - extract.flac",
%%	   "09 - David Kuckherman - Glen Velez Bendir technique - extract.flac",
%% 	   "10 - Stan Getz, João Gilberto - The Girl from Ipanema - extract.flac",
%songs(1).name;
%[songs(1).folder, stimname];
%[songPath,'vocoded',stimname];

for i = 1:length(songs),
    
    %% Load sound file
    stimname = songs(i).name;
    stim = [songs(i).folder,'/', stimname];
    %%[x, fs] = auload("file.flac")

    %% Change working directory for accessing sound files
    %refdir = 'vocoded/';
    savedir = 'music/vocoded/';
    refdir = '';
    
    audiopath = [stim];
    [sound, fs] = audioread(audiopath);

    tfs = 44100;
    tx = linspace(0, length(sound)/fs, length(sound));
    if (fs ~= tfs)
	sound = resample(sound, tx, tfs);
	ofs = fs;
	fs = tfs;
    end

    mblue = [0, 0.4470, 0.7410];
    %%waveform(sound, fs, mblue, 'en');

    % Configure the vocoding parameters
    p = struct();
    %p.analysis_filters  = filter_bands([100, 7000], nbbands, fs, 'greenwood', 1.5); % Order is multiplied by 4, so 1.5 gives 6
    % Test various filter orders in order to decrease overlap between frequency channels
    %p.analysis_filters  = filter_bands([100, HFlimit], nbbands, fs, 'greenwood', 3, 4); %Basal shift of +4mm along the cochlea Order is multiplied by 4, so 1.5 gives 6
    p.analysis_filters  = filter_bands([100, HFlimit], nbbands, fs, 'greenwood', 3, 0); %No Basal shift. Order is multiplied by 4, so 1.5 gives 6
    p.synthesis_filters = p.analysis_filters;
    
    % access to filterbank properties in p.analysis_filters.center, p.analysis_filters.upper, p.analysis_filters.lower, p.analysis_filters.order, p.analysis_filtes.filter_type.
    
    
    p.envelope = struct();
    %p.envelope.method = 'hilbert'; % fc is useless
    p.envelope.method = 'low-pass';
    p.envelope.rectify = 'half-wave';
    p.envelope.order = 2; % This produces a filter of order 4. Order 3 is not possible here.
    p.envelope.fc = 128;
    
    p.synth = struct();
    if strcmp(ctype, 'noise')
        p.synth.f0carrier = 'noise';
        p.synth.filter_before = false;
        p.synth.filter_after  = true;
    else
        p.synth.carrier = 'sin';
        p.synth.filter_before = false;
        p.synth.filter_after  = false;
    end

    % Generate example signals
    % TODO filter Carriers
    % Note that Carriers is not filtered, it is a wideband white noise
    [yOut, fsOut, p, ModC, Env, Bands, Carriers, levels] = vocode(sound, fs, p);


    % Save vocoded signal to file
    % yOut
    %plot(yOut)
    output = ['simulations_HFlim-', sprintf('%d', HFlimit), '_', stimname, '_', ctype, '_', sprintf('%02d', nbbands), 'bands_Vocoded.flac'];
    %outfile = [savedir, sprintf('%02i', i), '_', output]
    outfile = [savedir, sprintf('%02d', nbbands), 'bands/', output]

    status = mkdir(fullfile([savedir, sprintf('%02d', nbbands), 'bands/']));

    audiowrite(outfile, yOut, tfs);
    finalSig = yOut;
    
	% Save single vocoded bands to separate files
    if (vocbands == true),
					  % ModC.*levels
    %size(ModC)
    %plot(ModC(:,1).*levels(1))
	ModC = ModC.*levels;
	ModC = ModC * 0.6 / max(max(ModC));
	for j = 1:nbbands
        output = ['simulations_', stimname, '_', ctype, '_', sprintf('%02d', nbbands), 'bands_ModulatedCarrier_', sprintf('%02d', j), '_', sprintf('%02d', nbbands), '.flac'];
        outfile = [savedir, sprintf('%02d', nbbands), 'bands/', output]

        %outfile = [savedir, output]
        audiowrite(outfile, ModC(:,j), tfs);
	end
    end    

    % Save single original bands to separate files
    if (singleorig == true),
	%% Bands.*levels
	%% size(Bands)
	%% plot(Bands(:,1).*levels(1))
	Bands = Bands.*levels;
	Bands = Bands * 0.6 / max(max(Bands));
	for j = 1:nbbands
        output = ['simulations_', stimname, '_', ctype, '_', sprintf('%02d', nbbands), 'bands_OriginalBand_', sprintf('%02d', j), '_', sprintf('%02d', nbbands), '.flac'];
        outfile = [savedir, sprintf('%02d', nbbands), 'bands/', output]
        %outfile = [savedir, output]
        audiowrite(outfile, Bands(:,j), tfs);
	end
    end    
    
    %% Save recombined original bands to file
    if (recombinedorig == true),
	%% sum(Bands.*levels)
	%% plot(sum(Bands.*levels, 2))
	Bands = Bands.*levels;
	Bands = Bands * 0.6 / max(max(Bands));
	target = sum(Bands.*levels, 2);
	target = target * 0.6 / max(max(target));
	output = ['simulations_', stimname, '_', ctype, '_', sprintf('%02d', nbbands), 'bands_RecombinedOriginal', '.flac'];
		     %outfile = [savedir, sprintf('%02i', i), '_', output]
    outfile = [savedir, sprintf('%02d', nbbands), 'bands/', output]
    %outfile = [savedir, output]
	audiowrite(outfile, target, tfs);
    end

    %% Save carriers to separate files
    if (carriers == true),
	%% Carriers.*levels
	%% plot(Carriers(:,1).*levels(1))
	Carriers = Carriers.*levels;
	%% Carrier should be filtered to the same channel if ctype == "noise"
	%% TODO
	Carriers = Carriers * 0.6 / max(max(Carriers));
	for j = 1:nbbands
        output = ['simulations_', stimname, '_', ctype, '_', sprintf('%02d', nbbands), 'bands_Carrier_', sprintf('%02d', j), '_', sprintf('%02d', nbbands), '.flac'];
        outfile = [savedir, sprintf('%02d', nbbands), 'bands/', output]
        %outfile = [savedir, output]
        audiowrite(outfile, Carriers(:,j), tfs);
	end
    end
    
    %% Save envelopes to separate files
    if (envs == true),
	%% Env.*levels
	%% plot(Env(:,1).*levels(1));
	myEnv = Env.*levels;
	myEnv = 20*log10(myEnv./max(myEnv));
	%% plot(myEnv(:).*levels);
	myEnv = myEnv.*levels;
	myEnv = myEnv * 0.6 / max(max(myEnv));
	for j = 1:nbbands
        output = ['simulations_', stimname, '_', ctype, '_', sprintf('%02d', nbbands), 'bands_Envelope_', sprintf('%02d', j), '_', sprintf('%02d', nbbands), '.flac'];
        outfile = [savedir, sprintf('%02d', nbbands), 'bands/', output]
        %outfile = [savedir, output]
        audiowrite(outfile, myEnv(:,j), tfs);
	end
    end
    

end


