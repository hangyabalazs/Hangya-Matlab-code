function ibtwsubjcorr2_call
%IBTWSUBJCORR2_CALL    Calls IBTWSUBJCORR2.
%
%   See also IBTWSUBJCORR2.

pno = [{'31'} {'37'} {'39'} {'40'} {'n1'}];
pt = [{'virag'} {'lukacs'} {'gaal'} {'mojzsis'} {'wittner'}];
egs{1} = [{'177'} {'179'} {'181'}];
egs{2} = [{'12b'} {'12c'} {'24'}];
egs{3} = [{'10'} {'14'} {'40'}];
egs{4} = [{'125'} {'46'} {'46b'}];
egs{5} = [{'148a'} {'148b'} {'157'}];
for k1 = 1:4
    for k2 = k1+1:5
        ibtwsubjcorr2(pno{k1},pt{k1},egs{k1},pno{k2},pt{k2},egs{k2})
    end
end