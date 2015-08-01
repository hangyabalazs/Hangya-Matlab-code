ff = [inpdir1 bob '\' fname];       % load
load(ff)
unit = data(:,1)';

locunit = unit(ind1(k):ind2(k));
unit2 = locunit(1:const:end);

% hold on;plot(unit2,'r')

hold on;plot(linspace(1,length(unit2),length(locunit)),locunit,'c')