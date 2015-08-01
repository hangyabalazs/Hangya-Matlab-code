function wpcontrol
%WPCONTROL   Finger tapping test.
%   WPCONTROL performs a finger sequence tapping task to test procedural
%   memory. A five-element sequences are presented on the monitor, and the 
%   subject has to type the sequence using the non-dominant hand. Working
%   memory effects are excluded by presenting the sequence continuously. No
%   feedback on the typed characters is given to prevent interference.
%   The task consists of 30 sec. training periods interrupted by 30 sec.
%   breaks.
%
%   Reference: 
%   Marshall L, Helgadóttir H, Mölle M, Born J (2006) Boosting slow
%   oscillations during sleep potentiates memory. Nature 444:610-613.
%
%   See also WPTEST and WPRECALL2.

% Results directory
resdir = 'D:\MATLAB_R2007a\work\Balazs\Cognitive\wpres\';
cd(resdir)

% Time
ct = clock;
sct = [num2str(ct(1)) '_' num2str(ct(2)) '_' num2str(ct(3)) '_' ...
    num2str(ct(4)) '_' num2str(ct(5)) '_' num2str(round(ct(6)))];

% Screen size
sc = get(0,'ScreenSize');
scrn_width = sc(3);     % screen width
scrn_hight = sc(4);     % screen hight
switch scrn_width       % screen resolution code
    case 640
        scrn_res = 1;
    case 800
        scrn_res = 2;        
    case 1024
        scrn_res = 3;
    case 1152
        scrn_res = 4;
    case 1280
        scrn_res = 5;
    case 1600
        scrn_res = 6;
end

% Initialize Cogent
cogstd('spriority','high');     % increasing the priority class to 'high'
cgloadlib;                      % load all graphics libraries
config_display(1, 3, [0 0 0], [1 1 1], 'Arial', 32, 12); % configure display:
                                % full screen mode, resolution det. by monitor,
                                % black background, white foreground, 32 
                                % point Arial font, 12 offscreen buffers
config_results(['wpcontrol_' sct '.txt']); % configure result file
config_keyboard(100,5,'exclusive'); % configure keyboard: quelength = 100,
                                % resolution = 5, mode = 'exclusive'
start_cogent;                   % initialise Malab for running Cogent

% Create background
cgmakesprite(1,scrn_width,scrn_hight, 0.87, 0.92, 0.98);    % light blue background
cgsetsprite(1);       % destination for drawing commands: onscreen
cgpencol(0,0,0);     % font color
cgfont('Arial', 36); % font style and size
cgalign('c','c');    % align to center

% Task description
cgtext('A következõkben minden alkalommal öt karaktert fog látni a képernyõn.',0,250);
cgtext('Billentyûzze be a képernyõn látott karaktereket',0,150);
cgtext('a megadott sorrendben, a lehetõ leggyorsabban!',0,50);
cgtext('Csak a bal kezét használja!',0,-50);
cgtext('Szóljon, ha kezdhetjük!',0,-250);
cgsetsprite(0);       % destination for drawing commands: offscreen
cgdrawsprite(1,0,0);  % draw now
cgflip(0,0,0);        % copy the offscreen buffer to screen and then clear it
keyout = waitkeydown(inf, [52, 71]);
if keyout == 52    % escape
    stop_cogent;
    return
end

% Finger sequence tapping task
keylist = [1:26 28:36 63];
keylist5 = repmat(keylist,1,5);
lk = length(keylist);
snn = 2;    % number of periods
hitlists = cell(1,snn);
for sn = 1:snn   % training periods
    et = 0;
    hitlist = [];
    t1 = clock;
    while et < 30
        rp = randperm(lk*5);        % present characters
        cgmakesprite(2,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
        cgsetsprite(2);
        cgfont('Arial',60);
        tx = [keytable(keylist5(rp(1))) ' - ' keytable(keylist5(rp(2))) ' - '...
            keytable(keylist5(rp(3))) ' - ' keytable(keylist5(rp(4))) ' - '...
            keytable(keylist5(rp(5)))];
        cgtext(tx,0,0);
        cgsetsprite(0);
        cgdrawsprite(2,0,0);
        t = cgflip(0,0,0);

        keyout1 = waitkeydown(inf);     % key presses
        keyout2 = waitkeydown(inf);
        keyout3 = waitkeydown(inf);
        keyout4 = waitkeydown(inf);
        keyout5 = waitkeydown(inf);

        kin = [keylist5(rp(1)) keylist5(rp(2)) keylist5(rp(3)) keylist5(rp(4)) keylist5(rp(5))];
        kout = [keyout1 keyout2 keyout3 keyout4 keyout5];
        ishit = isequal(kin,kout);
        hitlist(end+1) = ishit;     % compare presented and pressed keys

        cgmakesprite(1,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
        cgdrawsprite(1,0,0);
        t = cgflip(0,0,0);
        pause(1)
        %     addresults('Number of recalled words:', hitno)
        et = etime(clock,t1);
    end
    cgmakesprite(1,scrn_width,scrn_hight, 0.87, 0.92, 0.98);    % light blue background
    cgsetsprite(1);       % destination for drawing commands: onscreen
    cgfont('Arial',36);
    if sn < snn
        cgtext('Fél perc után újabb teszt következik!',0,150);
    else
        cgtext('Köszönjük a közremûködést, vége a feladatnak!',0,150);  % end screen
    end
    cgsetsprite(0);       % destination for drawing commands: offscreen
    cgdrawsprite(1,0,0);  % draw now
    cgflip(0,0,0);        % copy the offscreen buffer to screen and then clear it
    if sn < snn
        pause(30)
    else
        pause(5)
    end
    
    hitlists{sn} = hitlist;
end

% Save recall lists
fnhl = ['hit_lists_' sct '.mat'];
save(fnhl,'hitlists');

% Stop Cogent
stop_cogent

% -------------------------------------------------------------------------
function str = keytable(K)
% Return strings for keys.

switch K
    case 1
        str = 'a';
    case 2
        str = 'b';
    case 3
        str = 'c';
    case 4
        str = 'd';
    case 5
        str = 'e';
    case 6
        str = 'f';
    case 7
        str = 'g';
    case 8
        str = 'h';
    case 9
        str = 'i';
    case 10
        str = 'j';
    case 11
        str = 'k';
    case 12
        str = 'l';
    case 13
        str = 'm';
    case 14
        str = 'n';
    case 15
        str = 'o';
    case 16
        str = 'p';
    case 17
        str = 'q';
    case 18
        str = 'r';
    case 19
        str = 's';
    case 20
        str = 't';
    case 21
        str = 'u';
    case 22
        str = 'v';
    case 23
        str = 'w';
    case 24
        str = 'x';
    case 25
        str = 'z';
    case 26
        str = 'y';
    case 27
        str = 'ö';
    case 28
        str = '1';
    case 29
        str = '2';
    case 30
        str = '3';
    case 31
        str = '4';
    case 32
        str = '5';
    case 33
        str = '6';
    case 34
        str = '7';
    case 35
        str = '8';
    case 36
        str = '9';
    case 53
        str = 'ü';
    case 54
        str = 'ó';
    case 57
        str = 'õ';
    case 58
        str = 'ú';
    case 61
        str = 'é';
    case 62
        str = 'á';
    case 63
        str = '0';
    case 65
        str = 'û';
        
    otherwise
        str = [];
        return
end