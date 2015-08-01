function wprecall
%WPRECALL   Prenight recall phase of WPTEST.
%   WPRECALL teaches the word pairs presented by WPTEST. Cue words are
%   shown and the subject has unlimited time to recall the appropriate
%   response word. After recall the correct pair appears on the screen. If
%   a minimum of 60% correct answers is not reached, cue words are
%   presented in a newly randomized order.
%
%   Reference: 
%   Marshall L, Helgadóttir H, Mölle M, Born J (2006) Boosting slow
%   oscillations during sleep potentiates memory. Nature 444:610-613.
%
%   See also WPTEST and WPRECALL2.

% Word pairs
wpdir = 'D:\MATLAB_R2007a\work\Balazs\Cognitive';
wpname = [wpdir '\word_pairs.txt'];
[wp1 wp2] = textread(wpname,'%s %s');

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
config_results(['wprecall_' sct '.txt']); % configure result file
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
cgtext('A következõkben a szópárok elsõ tagját fogja látni.',0,250);
cgtext('Billentyûzze be a megjelenõ szó párját, majd nyomja le az Entert!',0,150);
cgtext('Ha nem emlékszik, nyomja le a szóközt!',0,50);
cgtext('Szóljon, ha kezdhetjük!',0,-250);
cgsetsprite(0);       % destination for drawing commands: offscreen
cgdrawsprite(1,0,0);  % draw now
cgflip(0,0,0);        % copy the offscreen buffer to screen and then clear it
keyout = waitkeydown(inf, [52, 71]);
if keyout == 52    % escape
    stop_cogent;
    return
end

% Cued recall
hitno = 0;
wlim = 2;   % number of learned words to reach
while hitno < wlim
    rp = randperm(46);
    hitno = 0;
    for k = 1:4
        cgmakesprite(2,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
        cgsetsprite(2);
        cgfont('Arial',60);
        cgtext(wp1{rp(k)+4},-200,0);
        cgsetsprite(0);
        cgdrawsprite(2,0,0);
        t = cgflip(0,0,0);
        recalled_word = [];
        while 1
            keyout = waitkeydown(inf);
            if keyout == 52    % escape
                stop_cogent;
                return
            end
            if keyout == 71 | keyout == 59   % space, enter
                cgmakesprite(1,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
                cgdrawsprite(1,0,0);
                t = cgflip(0,0,0);
                break
            end
            if keyout == 55    % backspace
                if ~isempty(recalled_word)
                    recalled_word = recalled_word(1:end-1);
                    cgmakesprite(2,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
                    cgsetsprite(2);
                    cgtext(wp1{rp(k)+4},-200,0);
                    cgtext(recalled_word,200,0);
                    cgsetsprite(0);
                    cgdrawsprite(2,0,0);
                    t = cgflip(0,0,0);
                end
                continue
            end
            recalled_word = [recalled_word keytable(keyout)];
            cgmakesprite(2,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
            cgsetsprite(2);
            cgtext(wp1{rp(k)+4},-200,0);
            cgtext(recalled_word,200,0);
            cgsetsprite(0);
            cgdrawsprite(2,0,0);
            t = cgflip(0,0,0);
        end
        cgmakesprite(2,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
        cgsetsprite(2);
        cgtext(wp1{rp(k)+4},-200,0);
        recw2 = recalled_word;
        recw2(recw2=='i') = 'í';
        if isequal(recalled_word,wp2{rp(k)+4}) || isequal(recw2,wp2{rp(k)+4})
            cgpencol(0,1,0)
            hitno = hitno + 1;
        else
            cgpencol(1,0,0)
        end
        cgtext(wp2{rp(k)+4},200,0);
        cgpencol(0,0,0)
        cgsetsprite(0);
        cgdrawsprite(2,0,0);
        t = cgflip(0,0,0);
        pause(4.9)
        cgmakesprite(1,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
        cgdrawsprite(1,0,0);
        t = cgflip(0,0,0);
        pause(0.1)
    end
    addresults('Number of recalled words:', hitno)
    
    cgmakesprite(1,scrn_width,scrn_hight, 0.87, 0.92, 0.98);    % light blue background
    cgsetsprite(1);       % destination for drawing commands: onscreen
    cgfont('Arial',36);
    cgtext(['Helyesen felidézett szavak száma: ' num2str(hitno)],0,250);
    if hitno < wlim
        cgtext('A helyesen felidézett szavak száma nem érte el a 28-at.',0,150);
        cgtext('Újabb teszt következik!',0,50);
    else
        cgtext('A helyesen felidézett szavak száma elérte a 28-at.',0,150);
        cgtext('Köszönjük a közremûködést, vége a feladatnak!',0,50);
    end
    cgsetsprite(0);       % destination for drawing commands: offscreen
    cgdrawsprite(1,0,0);  % draw now
    cgflip(0,0,0);        % copy the offscreen buffer to screen and then clear it
    pause(10)
end

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
    case 65
        str = 'û';
        
    otherwise
        str = [];
        return
end