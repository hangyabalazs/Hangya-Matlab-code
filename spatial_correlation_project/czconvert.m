function czconvert(inpdir)
%CZCONVERT   Converts excel files to mat files.
%   CZCONVERT(INPDIR) converts 'allneurons' and 'position' xls files into
%   mat files based on the names in the headerline.

% Convert neurons
mm = pwd;
cd(inpdir)
[X,head] = xlsread([inpdir 'allneurons.xls']);

for k = 1:length(head)
    name = [head{k} '.mat'];
    data = X(:,k);
    data = data(~isnan(data));
    save(name,'data')
end

% Convert position data
[X,head] = xlsread([inpdir 'position.xls']);
for k = 1:length(head)
    name = [head{k} '.mat'];
    data = X(:,k);
    data = data(~isnan(data));
    save(name,'data')
end
cd(mm)