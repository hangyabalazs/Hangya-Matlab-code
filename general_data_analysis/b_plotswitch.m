function b_plotswitch
%PLOTSWITCH   Keypress function for THETASELECTOR_BETA plots.
%   Switch between various cutting values (see THETASELECTOR_BETA for more information).
%       q - switches up one plot
%       w - switches down one plot
%
%   See also THETASELECTOR_BETA.

inp = get(gcf,'CurrentCharacter');
switch inp
case 'q'
    swtch('u')
case 'w'
    swtch('d')
end

% -----------------------------------------------------------------------------------
function swtch(s)
fig = gcf;
plt = findobj(fig,'Type','line');
vis = get(plt,'Visible');
if isempty(find(strcmp(vis,'off')))
    set(plt(11:14),'Visible','off')
else
    ind1 = find(strcmp(vis,'on'));
    ind1 = ind1(find(ind1>10));
    switch s
    case 'u'
        ind2 = ind1 + 1;
    case 'd'
        ind2 = ind1 - 1;
    end
    r = rem(ind2,5);
    if r == 0
        r = 5;
    end
    ind2 = 10 + r;
    set(plt(ind1),'Visible','off')
    set(plt(ind2),'Visible','on')
end