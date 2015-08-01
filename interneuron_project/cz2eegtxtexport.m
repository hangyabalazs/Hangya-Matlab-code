%%

ff = 'X:\In_Vivo\balazs\_analysis\Czurko2\EEG\jh18s005_CSC01.txt';
[unit0 unit] = textread(ff,'%f %f');

ff = 'X:\In_Vivo\balazs\_analysis\Czurko2\EEG\temp.txt';
unit = textread(ff,'%f');

out = [unit0 unit];
fo = 'X:\In_Vivo\balazs\_analysis\Czurko2\EEG\jh18s005_CSC12.txt';
save(fo,'out','-ASCII')