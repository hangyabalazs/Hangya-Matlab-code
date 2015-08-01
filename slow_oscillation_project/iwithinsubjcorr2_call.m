function iwithinsubjcorr2_call
%IWITHINSUBJCORR2_CALL    Calls IWITHINSUBJCORR2.
%
%   See also IWITHINSUBJCORR2.

pno = [{'31'} {'37'} {'39'} {'40'} {'n1'} {'n2'}];
pt = [{'virag'} {'lukacs'} {'gaal'} {'mojzsis'} {'wittner'} {'bak'}];
egs{1} = [{'177'} {'179'} {'181'} {'232'} {'285'}];
egs{2} = [{'12b'} {'12c'} {'24'} {'52'} {'53'}];
egs{3} = [{'10'} {'14'} {'15'} {'18'} {'40'}];
egs{4} = [{'125'} {'46'} {'46b'}];
egs{5} = [{'148a'} {'148b'} {'157'} {'262'} {'298'}];
egs{6} = [{'90a'} {'90b'} {'90c'}];
for k = 1:6
    iconvdivcorr3(pno{k},pt{k},egs{k})
end