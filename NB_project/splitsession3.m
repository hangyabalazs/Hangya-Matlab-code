%% Names and IDs

% What is the ID of the session you would like to split?
session = '150108a';

% Which tetrode do you wish to re-cluster?
tetrode = 'TT6';

%% Split point

% A split point has to be defined; a split point can be obtained by
% visualizing the tetrode data in MClust, choosing a split point and
% multiplying the MClust timestamp by 100 to be compatible with Neuralynx
% timestamps
ct = 1.9855 * 10^12;    % MClust timestamp * 100

% Do you need the segment before or after the split point?
ba = 'after';   % options: 'before' or 'after'
rba = '<';
if isequal(ba,'after')
    rba = '>';
end

%% Folders

% You need a folder for temporary storage ('junk folder'):
global DATAPATH
junkf = [DATAPATH 'NB\junk\'];

% Copy the session folder from your CellBase to your junk folder to split
% the session while keeping the original data in your CellBase!

% Rename the copied session by putting a ' - Copy' tag to avoid confusion
% of your folders!

% The name of your temporary folder is now this:
junkf2 = [junkf session ' - Copy\'];

% The full path of your raw tetrode recording is this:
fn = [junkf2 tetrode '.ntt'];


%% Read data

% Import Neuralynx tetrode data
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatSpike(fn,[1 1 1 1 1],1,1,1);

%% Restrict data to the desired segment

% We overwrite all imported variables with the restricted versions
eval(['inx = find(TimeStamp' rba 'ct);']);   % indices for the required time stamps
Samples = Samples(:,:,inx);
TimeStamp = TimeStamp(inx);
SampleFrequency = SampleFrequency(inx);
ChanNum = ChanNum(inx);
NumValSamples = NumValSamples(:,inx);

%% Write data

% We write the data to the same folder with an underscore after the tetrode name 
fn = [junkf2 tetrode '_.ntt'];   % file name
Mat2NlxTT(fn,0,1,1,length(TimeStamp),[1 1 1 1 1 1],...
    TimeStamp,ChanNum,SampleFrequency,NumValSamples,Samples,NlxHeader)   % export data as .ntt file

%% TrialEvents

% We also need to restrict the behavioral time stamps in the trial events file 
fn = [junkf2 'TrialEvents.mat'];  % file name
TE = load(fn);   % trial events structure
ctt = ct / 1000000;   % convert time stamp
eval(['tinx = find(TE.TrialStart' rba 'ctt);']);   % index set for the required time stamps

fld = fieldnames(TE);   % all fields has to be modified
for k = 1:length(fld)   % loop through the fields of the trial events structure
    TE.(fld{k}) = TE.(fld{k})(tinx);   % restrict the field
end

save(fn,'-struct','TE')   % overwrite trial events structure

%% Follow the instructions!

% Now delete all unnecessary files from your session in the junk folder!
% This includes 
%   - raw data except your newly created ntt file
%   - MClust feature files
%   - old cluster files
%   - old STIMSPIKES and EVENTSPIKES files

% Copy your session from the junk folder back to CellBase!

% Rename it with a single letter to make it a valid session code! (e.g. YYMMDDx)!

% Rename your ntt file - remove the underscore!

% Re-cluster your new ntt file in MClust!

% CAREFUL! quickanalysis overwrites TrialEvents!