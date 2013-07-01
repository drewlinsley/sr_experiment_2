function got_song(command,song)

if ~exist('song','var'),
   song = []; 
end
nrchannels = 1; % One channel only -> Mono sound.
freq = 8192;      % Fs is the correct playback frequency for handel.
switch command
    case 'initialize'
        InitializePsychSound;
    case 'play'
        %scene_sound = [A_note(whole_note),D_note(whole_note),E_note(half_note),F_note(half_note)];
        
        %GOT theme song:
        %scene_sound = [A_note(whole_note),D_note(whole_note),E_note(half_note),F_note(half_note),...
        %A_note(whole_note),D_note(whole_note),E_note(half_note),F_note(half_note),...
        %   A_note(whole_note),A_note(whole_note),A_note(whole_note),...
        %  D_note(whole_note),D_note(whole_note),D_note(whole_note),...
        % D_note(whole_note),E_note(half_note),F_note(whole_note),F_note(whole_note),...
        %D_note(whole_note),D_note(whole_note),D_note(whole_note)];
        pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
        PsychPortAudio('FillBuffer', pahandle, song);
        PsychPortAudio('Start', pahandle, [], 0, 1);
    case 'close'
        pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
        PsychPortAudio('Stop', pahandle);
        PsychPortAudio('Close', pahandle);
end



