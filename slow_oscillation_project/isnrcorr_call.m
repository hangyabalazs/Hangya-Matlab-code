function isnrcorr_call
%ISNRCORR_CALL    Calls ISNRCORR.
%
%   See also ISNRCORR.

pno = [{'31'} {'37'} {'39'} {'40'} {'n1'} {'n2'}];
pt = [{'virag'} {'lukacs'} {'gaal'} {'mojzsis'} {'wittner'} {'bak'}];
egs{1} = [{'177'} {'179'} {'181'} {'232'} {'285'}];
egs{2} = [{'12b'} {'12c'} {'24'} {'52'} {'53'}];
egs{3} = [{'10'} {'14'} {'40'} {'15'} {'18'}];
egs{4} = [{'125'} {'46'} {'46b'}];
egs{5} = [{'148a'} {'148b'} {'157'} {'262'} {'298'}];
egs{6} = [{'90a'} {'90b'} {'90c'}];
for k = 1:6
    for t = 1:length(egs{k})
        isnrcorr(pno{k},pt{k},egs{k}{t})
    end
end