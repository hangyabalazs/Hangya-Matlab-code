function [Xv, Xd, V, d]=lda(X,CLASSES)
%
% LDA   Linear Discriminant Analysis
%
%   [Xd, V, d]=lda(X,CLASSES)
%
%   IN: 
%       X is an NxM matrix with N samples of M element vectors
%       CLASSES contains for each class the pointers to the relevant samples.
%               either an IxJ matrix with I*J = N where I classes of J
%               samples each
%               or a cellaray of vectors with each vector representing a
%               class with the elements indexing X
%
%   OUT:
%       Xd is the transformation of X into the subpace found
%       V is the vector of discriminants
%       d is the % of variance accounted for by each dir
%

%% AK 2002

percent = 0.99999999999999999;

[nSamples,nVec] = size(X);
if ~iscell(CLASSES)
    nClasses = size(CLASSES,1);
    for iC = 1:nClasses
        class = CLASSES(iC,~isnan(CLASSES(iC,:)));
        C{iC} = class;
    end
    CLASSES = C;
    clear C;
end
nClasses = length(CLASSES);
sprintf('LDA: %d element vector, %d classes, %d samples \n',nVec,nClasses,nSamples)
%--------------------------
%
%Sw = zeros(nVec,nVec); 
%Sb = zeros(nVec,nVec);
%m_all=mean(X);
meansClass = zeros(nVec,nVec);
meansAll   = zeros(nVec,nVec);
classV = zeros(1,nSamples);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate within (Sw) and between (Sb) scatter matrices 
for (iClass=1:nClasses)
   class  = CLASSES{iClass};
   classV(class) = iClass; 
   meansClass(iClass,:)=mean(X(class,:));
end
meansAll = meansClass(classV,:);

Sb = cov(meansClass);

Sw = cov(X - meansAll);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the trick, generalized eigenvalue decompositon
%[V,D] = eig(Sw,Sb);
[V,D] = eig(Sw \ Sb);

%[U,sm,X,V,W] = cgsvd(A,L)

%[U,V,X,C,S] = gsvd(Sw,Sb);
%d=sqrt(diag(C'*C)./diag(S'*S));
%D=diag(d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now find out how many components are needed to account for most data
[d,ind]=sort(diag(D));    
ind = ind(~isinf(d)); d = d(~isinf(d));     % make sure sure inf's dont mock things up 
drank=abs(d)/sum(abs(d))*100;

%reverse
d = d(end:-1:1); 
ind = ind(end:-1:1); 
drank = drank(end:-1:1);  



cumdrank = cumsum(drank);
p = find(cumdrank > percent );
if isempty(p)
    sprintf('LDA: Only %.2f % variance is accounted for.',cumdrank(end))
    rank = length(cumdrank);
else
    rank = p(1);    % how many dims are needed to account for 'percent' of the data
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the important (lambda big) discriminant directions
V = V(:,ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now convert the data into the new coordinates
for (iV = 1:rank)
   v = V(:,iV);
   Xd(:,iV) = X*v;
end

for (iClass=1:nClasses)
   class = CLASSES{iClass};
   Xv{iClass,1} = Xd(class,1);
   if size(Xd,2) > 1
        Xv{iClass,2} = Xd(class,2);
   end
end