function b_shortan
%SHORTAN    Shortan filenames through leaving the last character.

cd d:\_analysis\matlab_data\Wavelet\segments_long
d = dir(pwd);
for k = 3:length(d)
    fname = d(k).name;
    ext = fname(end-3:end);
    fnamenew = [fname(1:end-5) ext];
    load(fname)
    eval(['save(fnamenew,''ThetaSegments'')']);
end