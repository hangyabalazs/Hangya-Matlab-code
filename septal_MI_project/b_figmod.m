function b_figmod
%FIGMOD Standardizes the IVSI figs.
%   FIGMOD sets the axes properties and the markers of the IVSI figs.
%
%   See also IVSI and IVSIRUN.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
where=input('Where are the files? ','s');
cd(where);
files=dir(where);
files=files(3:end);
sf=length(files);
for i=1:sf
   fnm=files(i).name;
   ffnm=[where fnm];
   open(ffnm)

% Standardize figure properties   
   H = findobj(gcf,'Color','r','Marker','o');
   set(H,'Marker','s');
   axis([0 450 0 800]);
   line([0 450],[0 450],'Color','g');
   saveas(gcf,fnm,'fig');
   close
end;