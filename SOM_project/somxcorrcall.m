function somxcorrcall(ccData)
%SOMXCORRCALL   Caller function for SOMXCORR.
%   SOMXCORRCALL calls SOMXCORR for cross-correlatin calculation on a
%   sequence of files.
%
%   See also SOMXCORR.

% Calling 'somxcorr'
for k = 15:length(ccData)
    spkd = ccData{k}(1).SpkData;
    try
        somxcorr(spkd);
    catch
        disp(num2str(k))
        lasterr
    end
end