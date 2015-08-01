function chirpvec = chordchirp(freq_vec,attn_vec,sampf,chirp_dur,dirn,sweep_speed,sweeptype)

 

if nargin < 1,

    freq_vec = [4000:2000:20000];

end

if nargin < 2,

    attn_vec = 1*ones(length(freq_vec),1);

end

if nargin < 3,

    sampf = 40000;

end

if nargin < 4,

    chirp_dur = 0.5;

end

if nargin < 5,

    dirn = 'up';

end

if nargin<6,

    sweep_speed = 5000; % Hz/s

end

if nargin<7,

    sweeptype = 'quadratic';

end

 

RaiseFallDuration = 0.08;

minfreq = 1000;

maxfreq = 20000;

switch dirn

    case 'up'

        f1 = [freq_vec + (sweep_speed/chirp_dur)];

        f1(f1>maxfreq) = maxfreq;

    case 'down'

        f1 = [freq_vec-(sweep_speed/chirp_dur)];

        f1(f1<minfreq) = minfreq;

end

 

% Make the other freq vec.

timevec = (0:1/sampf:chirp_dur)';

chirpmat = nan(length(freq_vec),length(timevec));

 

for iF = 1:length(freq_vec),

    switch dirn

        case 'up'

            chirpmat(iF,:) = chirp(timevec,freq_vec(iF),chirp_dur,f1(iF),sweeptype).*attn_vec(iF);

        case 'down'

            chirpmat(iF,:) = chirp(timevec,freq_vec(iF),chirp_dur,f1(iF),sweeptype).*attn_vec(iF);

    end

end

 

chirpvec = sum(chirpmat,1);

 

chirpvec = apply_risefall(chirpvec,RaiseFallDuration,sampf);

 

function SoundWaveform = apply_risefall(SoundWaveform,RaiseFallDuration,SamplingRate)

 

TimeVec = (0:1/SamplingRate:RaiseFallDuration)';

RaiseVec = linspace(0,1,length(TimeVec));

 

if(length(RaiseVec)<length(SoundWaveform))

    SoundWaveform(1:length(TimeVec)) = RaiseVec.*SoundWaveform(1:length(TimeVec));

    SoundWaveform(end-length(TimeVec)+1:end) = RaiseVec(end:-1:1).*SoundWaveform(end-length(TimeVec)+1:end);

else

    warning('Sound length is too short to apply rise and fall envelope');

end

return