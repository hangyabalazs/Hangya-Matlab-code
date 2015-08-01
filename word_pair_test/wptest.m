function wptest
%WPTEST   Word pair test.
%   WPTEST performs the learning phase of a word pair test (declarative 
%   memory test). Semantically related pairs of nouns (46) are presented in
%   a randomized order (1/5 s presentation rate, 100 ms interstim. iv.). 
%   Four dummy pairs at the beginning and end are presened to buffer
%   primacy and recency effects.
%
%   Reference: 
%   Marshall L, Helgadóttir H, Mölle M, Born J (2006) Boosting slow
%   oscillations during sleep potentiates memory. Nature 444:610-613.
%
%   See also WPRECALL and WPRECALL2.

% Word pairs
wpdir = [matlabroot '\work\Balazs\word_pair_test'];
wpname = [wpdir '\word_pairs.txt'];
[wp1 wp2] = textread(wpname,'%s %s');

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
% config_results(['wptest_' sct '.txt']); % configure result file
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
cgtext('A következõkben 54 szópárt fog látni.',0,250);
cgtext('Próbálja meg megjegyezni az összetartozó szavakat!',0,150);
cgtext('Szóljon, ha kezdhetjük!',0,-250);
cgsetsprite(0);       % destination for drawing commands: offscreen
cgdrawsprite(1,0,0);  % draw now
cgflip(0,0,0);        % copy the offscreen buffer to screen and then clear it
keyout = waitkeydown(inf, [52, 71]);
if keyout == 52    % escape
    stop_cogent;
    return
end

% Dummy I.
rp = randperm(4);
for k = 1:4
    cgmakesprite(2,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
    cgsetsprite(2);
    cgfont('Arial', 60);
    cgtext(wp1{rp(k)},-200,0);
    cgtext(wp2{rp(k)},200,0);
    cgsetsprite(0);
    cgdrawsprite(2,0,0);
    t = cgflip(0,0,0);
    pause(4.9)
    
    cgmakesprite(1,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
    cgdrawsprite(1,0,0);
    t = cgflip(0,0,0);
    pause(0.1)
end

% 46 pairs
rp = randperm(46);
for k = 1:4
    cgmakesprite(2,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
    cgsetsprite(2);
    cgtext(wp1{rp(k)+4},-200,0);
    cgtext(wp2{rp(k)+4},200,0);
    cgsetsprite(0);
    cgdrawsprite(2,0,0);
    t = cgflip(0,0,0);
    pause(4.9)
    
    cgmakesprite(1,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
    cgdrawsprite(1,0,0);
    t = cgflip(0,0,0);
    pause(0.1)
end

% Dummy II.
rp = randperm(4);
for k = 1:4
    cgmakesprite(2,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
    cgsetsprite(2);
    cgtext(wp1{rp(k)+50},-200,0);
    cgtext(wp2{rp(k)+50},200,0);
    cgsetsprite(0);
    cgdrawsprite(2,0,0);
    t = cgflip(0,0,0);
    pause(4.9)
    
    cgmakesprite(1,scrn_width,scrn_hight, 0.87, 0.92, 0.98);
    cgdrawsprite(1,0,0);
    t = cgflip(0,0,0);
    pause(0.1)
end

% End screen
cgmakesprite(1,scrn_width,scrn_hight, 0.87, 0.92, 0.98);    % light blue background
cgsetsprite(1);       % destination for drawing commands: onscreen
cgfont('Arial', 36);
cgtext('Véget ért a szavak bemutatása.',0,250);
cgtext('Köszönjük a figyelmet!',0,50);
cgsetsprite(0);       % destination for drawing commands: offscreen
cgdrawsprite(1,0,0);  % draw now
cgflip(0,0,0);        % copy the offscreen buffer to screen and then clear it
pause(10)

% Stop Cogent
stop_cogent