function chirpvec = sweep(freq1,freq2,sampf,chirp_dur)

% Time vector
t = (0:1/sampf:chirp_dur);

% Chirp
y = chirp(t,freq1,chirp_dur,freq2);

% Apply rise and fall attenuation
RaiseFallDuration = 0.08;
chirpvec = apply_risefall(y,RaiseFallDuration,sampf);

% -------------------------------------------------------------------------
function SoundWaveform = apply_risefall(SoundWaveform,RaiseFallDuration,SamplingRate)

TimeVec = (0:1/SamplingRate:RaiseFallDuration)';
RaiseVec = linspace(0,1,length(TimeVec));
if length(RaiseVec) < length(SoundWaveform)
    SoundWaveform(1:length(TimeVec)) = RaiseVec .* SoundWaveform(1:length(TimeVec));
    SoundWaveform(end-length(TimeVec)+1:end) = RaiseVec(end:-1:1) .* ...
        SoundWaveform(end-length(TimeVec)+1:end);
else
    warning('Sound length is too short to apply rise and fall envelope');
end