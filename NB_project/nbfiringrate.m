%% all ChAT and pChAT cells

ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified ChAT+ cells
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative ChAT+ cells
allChAT = [ChAT pChAT];

%% n029_120202a_3.5
% Firing rate is very stable despite variable behaviour. There is not much
% data outside the task; the tracking look useless; it was probably all
% recorded in headfixed, starts with the task in 12 s after recording
% start, stimulation sessions are after the task. ~1 min after the stim
% session; postsim. FR is slightly LOWER before end of session.

firingrate_analysis(allChAT{2})

%% n029_120203a_3.1
% The firing rate in the session is steadily DECREASING from 6 to 3 (the
% cluster is good and not drifting into noise, in fact, this is used as
% example cluster for identified cholinergic neurons). The cell is similar
% to the previous one (but it is another cell); the behaviour is also
% similar, from impulsive through good to unmotivated. Based on tracking,
% the whole session is headfixed. It starts with the session and tagging
% comes right after. In the ~1 min poststim. epoch (still fixed), the FR is
% HIGHER.

firingrate_analysis(allChAT{3})

% SLEEP
% Firing seem to be correlated w sleep: almost no firing in sleep (no REM
% recorded). LFP is very noisy during movement; delta-to-all-power ratio
% shows sleep. FR does not seem to correlate with speed beyond what is
% explained by sleep.
sleepsess = getvalue('alternative_sessions',allChAT{3});
firingrate_analysis_sleep(sleepsess{1})


%% n029_120214b_2.2
% Initially lower firing rate is probably a clustering issue. Pre-session
% segment only 9s, not useful. Post-stim. FR seems LOWER (only 25s). No
% other FR change. (Session structure: behav. session, ~45s gap,
% stimulation - prob. headfixed, 25s before end of rec. Behav. looks
% stable.)

firingrate_analysis(allChAT{4})

% The same cell in the next day session
% Session seem to be all headfixed. Slow change in FR may be explained by
% cluster drifting. FR does not change post-session or post-stim (both
% ~1min). Session structure is the same as for the prev. session. Behav. is
% somewhat less stable, imp. in the beginning.
nextsess = getvalue('alternative_sessions',allChAT{4});
firingrate_analysis(nextsess{1})

%% n029_120301a_3.3
% FR DECREASE throughout the session is not explaind by cluster drift: the
% cluster drifts, but in the opposite way. Instead, it seems to correlate
% with the hit rate. Same session structure: only headfixed, start w
% behaviour, then tagging. Post-session segment (~1min) does not show a
% change beyond the overall decrease throughout the session. Post-stim FR
% (23s) is HIGHER.

firingrate_analysis(allChAT{5})

%% n046_121229a_1.1
% Stable cluster. The animal is not too motivated for the first 100 trials,
% then suddenly becomes impulsive. At that transition, the FR also
% INCREASES (even for miss/CR trials) when calculated for full trials.
% Before the tone, it rather decreases if anything.
% The session starts w freely moving tagging. FR is much HIGHER in freely
% moving: 5.6 Hz both pre-stim. (16s) and post-stim (55s). No segments
% after the behavior (next behav. session starts in the next recording
% right after).

firingrate_analysis(allChAT{6})

% Behaviour session without background noise
% Continuation of prev. session, but w/o noise. All headfixed. FR is LOWER 
% post-session (46s). (Post.-stim also, but only 6s.)
nextsess = getvalue('alternative_sessions',allChAT{6});
firingrate_analysis(nextsess{1}{1})

% Behaviour session without feedback-delay
% Part of the last session (same day). No stimulation file. FR is
% DECREASING slowly; not due to drift (checked in MClust as well).
firingrate_analysis(nextsess{1}{2})


%% n046_121230a_1.2
% Firing rate is steadily DECREASING; probably not due to drift (MClust);
% maybe correlated w hit rate? (both hit and FA rate decreases slowly).
% Session structure: freely moving tagging, behav., tagging again
% (headfixed?). FR is substantially HIGHER freely moving (both before and
% after tagging).

firingrate_analysis(allChAT{7})

%% n046_121231a_1.2
% Slow FR change may be because of cluster drift. Freely moving tagging
% before the session: pre-stim. FR somewhat higher, post-stim. FR 
% definitely HIGHER. After session, head-fixed tagging; not much FR change
% there.

firingrate_analysis(allChAT{8})

%% n046_130101a_6.1
% At around ~50 trials, the hit rate dips and FA rate peaks; it is
% accompanied by a FR INCREASE. FR is LOWER after the animal stops. Session
% starts w freely moving tagging. FR is much HIGHER in freely moving, both
% pre- and post-stim. (17 and 73s). The session ends with head-fixed
% tagging. FR somewhat HIGHER post-stim., but lower then in freely moving.

firingrate_analysis(allChAT{9})

% No noise session
% All headfixed. Tagging after behavior. FR somewhat HIGHER post-stim.
% (same as the first session).
nextsess = getvalue('alternative_sessions',allChAT{9});
firingrate_analysis(nextsess{1}{1})

% No feedback-delay session
% FR is somewhat HIGHER and MORE VARIABLE when the animal is not doing the
% task (have to check whether that's because of the shorter trials). All
% headfixed; tagging before and after behav.
firingrate_analysis(nextsess{1}{2})

% SLEEP
firingrate_analysis_sleep(nextsess{1}{3})

% Next day session
% A dip in hit rate now not accompanied by FR change; stable FR throughout
% the session. Same: freely moving FR much HIGHER; poststim. FR somewhat
% HIGHER after headfixed tagging. (Session: freely moving tagging, behav.,
% headfixed tagging.)
firingrate_analysis(nextsess{1}{4})

% No noise session, not doing
% All headfixed. Fairly HIGH FR for this cell.
firingrate_analysis(nextsess{1}{5})

% SLEEP
firingrate_analysis_sleep(nextsess{1}{6})

% No-FD session
% Stable (relatively low) FR throughout the session. Same: freely moving FR
% much HIGHER; poststim. FR somewhat HIGHER after headfixed tagging.
% (Session: freely moving tagging, behav., headfixed tagging.)
firingrate_analysis(nextsess{1}{7})

% No noise session
% All headfixed. FR does not change w behav. (hit of FA rate). FR HIGHER
% after first headfixed tagging; not after the second (only 10s).
firingrate_analysis(nextsess{1}{8})

% Next day session
% Session: freely moving tagging, behav., headfixed tagging. FR is somewhat
% decreasing around the end of the behav., can be clustering issue, but FR
% comes back. FR is HIGHER on freely moving and also after headfixed
% tagging.
firingrate_analysis(nextsess{1}{9})

%% n046_130102x_4.1
% Session: freely moving tagging, behav., headfixed tagging; stable FR. We
% don't have the cell for the freely moving part. It does not change FR
% after the animal stops or after tagging.

firingrate_analysis(allChAT{10})

% No noise; not doing
% No changes in FR.
nextsess = getvalue('alternative_sessions',allChAT{10});
firingrate_analysis(nextsess{1}{1})

%% n046_130103a_6.2
% Session: freely moving tagging, behav., headfixed tagging. FR DECREASES
% throughtout the session (small cell, but prob. not clustering issue based
% on MClust). However, FR gets much HIGHER after the animal stops
% performing. FR is also HIGHER in freely moving.

firingrate_analysis(allChAT{11})

%% n046_130104a_6.2
% Session: freely moving tagging, behav., headfixed tagging. FR DECREASES
% throughtout the seesion (small cell, hard to say whether it is clustering
% issue from on MClust). FR HIGHER in freely moving as compared to behav.
% session average, but not much higher than the beginning of the behavioral
% session.

firingrate_analysis(allChAT{12})

% SLEEP
% FR is only slightly LOWER in sleep. It does not seem to be correlated w
% sleep either. However, FR has a sharp peak when the animal first wakes
% up.
nextsess = getvalue('alternative_sessions',allChAT{12});
firingrate_analysis_sleep(nextsess{1}{1})

%% n046_130108a_4.1
% Session: behav., tagging; all headfixed. FR is slightly INCREASING
% throughout the session (prob. not drifting based on MClust). FR is
% slightly HIGHER poststim.

firingrate_analysis(allChAT{13})

% SLEEP
% FR somewhat LOWER in sleep and somewhat correlated w speed; not too
% strong.
nextsess = getvalue('alternative_sessions',allChAT{13});
firingrate_analysis_sleep(nextsess{1}{1})

% No-FD session
% Session: freely moving tagging, behav., headfixed tagging. FR is
% INCREASING throughout the session (prob. not drift based on MClust). FR
% keeps INCREASING after the session. FR is not really diff. in freely
% moving.
firingrate_analysis(nextsess{1}{2})

%% n046_130108a_4.2
% Session: behav., tagging; all headfixed. Very low FR (nice cluster, no
% drift based on MClust). 

firingrate_analysis(allChAT{14})

% SLEEP
% FR is somewhat HIGHER in freely moving (still <1Hz).
nextsess = getvalue('alternative_sessions',allChAT{14});
firingrate_analysis_sleep(nextsess{1}{1})

% SLEEP #2
% FR does not change (awake seg. is very short though); the higher FR
% segment correspond to stimulation (stim. does not do anything to the LFP,
% even though this is in sleep).
firingrate_analysis_sleep(nextsess{1}{2})

% No noise session
% All headfixed; nothing interesting.
firingrate_analysis(nextsess{1}{3})

% GO tone changed
% Cell (and another) resonds to go tone the same strong way.
firingrate_analysis(nextsess{1}{4})

% No-FD session
% Session: freely moving tagging, behav., headfixed tagging. Maybe HIGHER
% FR after 2nd tagging. FR is not diff. for freely moving. 
firingrate_analysis(nextsess{1}{5})

%% 'n045_121217x_4.6'
% Low FR. Seems to be correlated with hit and FA rate.

firingrate_analysis('n045_121217x_4.6')

%% pChATs - n018_111018a_7.1
% All headfixed. Sudden increase in FR does not corr. w behavior; no
% clustering problem either. No interesting FR change.

firingrate_analysis(allChAT{15})

%% pChATs - n023_111220a_1.2
% All headfixed. Strong DECREASE of FR throughout the session; can be in
% part clustering, but unlikely that it is only drift. Poststim FR (after
% behav.) is somewhat HIGHER (note that this is a pChAT cell), nut the seg.
% is only 15s.

firingrate_analysis(allChAT{16})

%% pChATs - n029_120207b_1.1
% All headfixed. FR DECREASING in session.

firingrate_analysis(allChAT{17})

%% pChATs - n029_120210a_3.3
% All headfixed. FR slowly DECREASES but goes to NEAR 0 when the animal
% stops (at a point; hard to say when does the animal 'gives up'); perfect
% cluster. maybe slightly HIGHER post-stim (after behav)?

firingrate_analysis(allChAT{18})

%% pChATs - n028_120211a_8.1
% All headfixed. FR correlated w engagement: FR slightly INCREASES when the
% animal starts; then DECREASES throughout the session. Stop segment: FR is
% NEAR 0; but post-session (1 min) HIGH again.

firingrate_analysis(allChAT{19})

%% pChATs - n028_120211a_8.2
% All headfixed. EXACT SAME PATTERN as prev. cell (same session!!!).

firingrate_analysis(allChAT{20})

%% pChATs - n029_120215a_3.4
% All headfixed. No FR changes.

firingrate_analysis(allChAT{21})

%% pChATs - n029_120220a_3.1
% All headfixed. Atypically high FR. No changes, except stop seg is LOW
% (only 15s) and post-stim is LOW.

firingrate_analysis(allChAT{22})

% Second session right after the first
% Relatively stable, high FR.
nextsess = getvalue('alternative_sessions',allChAT{22});
firingrate_analysis(nextsess{1}{1})

%% pChATs - n029_120221b_6.1
% All headfixed. FR maybe DECREASING; somewhat corr. w. FA rate.
% Postsession slightly LOWER?

firingrate_analysis(allChAT{23})

%% pChATs - n029_120222b_4.1
% All headfixed. FR increase then decrease; does not seem to corr. w.
% behav. too much. Postsession somewhat HIGHER.

firingrate_analysis(allChAT{24})

%% pChATs - n029_120313a_1.1
% All headfixed. FR is correlated w engagement/impulsivity/FA rate/movement
% (INCREASE, then DECREASE). Pre-session LOW; post-session somewhat HIGHER.

firingrate_analysis(allChAT{25})

%% pChATs - n029_120314a_3.1
% All headfixed. FR changes not corr w. behav.; maybe drift - small cell.

firingrate_analysis(allChAT{26})

%% pChATs - n037_121006a_4.1
% All headfixed. FR is corr w hit rate: INCREASE then DECREASE.
% Post-session somewhat HIGHER; post-stim. (after session) even HIGHER.

firingrate_analysis(allChAT{27})

%% pChATs - n046_121210a_8.1
% All headfixed. Sudden increase in hit rate is accompanied by an INCREASE
% in FR; then slow DECREASE; sudden decrease of hit rate before stop is not
% reflected in FR.

firingrate_analysis(allChAT{28})

%% pChATs - n046_121213a_3.1
% Session: freely moving tagging, behav., headfixed tagging. First postim.
% FR very HIGH (freely moving), second not. Unfortunately, freely moving
% prestim. segment is very short: 5s (FR normal there). FR is DECREASING:
% it is higher while the animal is imp. (first half) and lower when
% unmotivated (second half); not so clear from pre-tone segments though.

firingrate_analysis(allChAT{29})

%% pChATs - n046_121218a_2.2
% All headfixed. FR might be DECREASING: higher when impu;sive, lower when
% not doing; not so clear from pre-tone segments though. FR LOWER when
% post-session.

firingrate_analysis(allChAT{30})

%% pChATs - n046_121219a_8.1
% All headfixed. FR might be DECREASING: higher when impu;sive, lower when
% not doing; not so clear from pre-tone segments though. FR LOWER
% post-session and in stop segment; gets somewhat HIGHER again post-stim.

firingrate_analysis(allChAT{31})

%% pChATs - n045_121231a_8.1
% All headfixed; FR is DECREASING, maybe paralelling behav. (but not clear
% that initial not-doing part is accompanied by lower FR). Stop segment
% maybe slightly LOWER FR.

firingrate_analysis(allChAT{32})

%% pChATs - n046_130102x_4.3
% THIS IS PART OF THE SESSION: POST-STIM1 IS COMPROMISED. FR is INCREASING,
% reflecting increasing involvement!!!

firingrate_analysis(allChAT{33})

% No noise session
% All headfixed. Short; FR may reflect hit rate; FA rate is constant zero.
nextsess = getvalue('alternative_sessions',allChAT{33});
firingrate_analysis(nextsess{1}{1})

%% pChATs - n046_130104a_4.1
% FR increasing due to cluster drift.

firingrate_analysis(allChAT{34})

% 3 days later
% FR increasing, but can be cluster issue, small cell.
nextsess = getvalue('alternative_sessions',allChAT{34});
firingrate_analysis(nextsess{1}{1})

%% pChATs - n046_130104a_6.1
% Session: freely moving tagging, behav., headfixed tagging. Does not seesm
% to fire much more in freely moving (not moving much though). FR INCREASES
% in the second half of the behav. session; no apparent chenge in behavior.

firingrate_analysis(allChAT{35})

%% pChATs - n046_130108a_8.1
% All headfixed. FR is DECREASING with decreaing involvement, but once the
% mouse decided to stop, it INCREASES.

firingrate_analysis(allChAT{36})

% Previous day
% All headfixed. This time, FR is quite stable (maybe slight decrease?);
% again increase after the mouse stops (short seg.0 - although the number
% shows different, because it is only true for the second half of the stop
% segment.
nextsess = getvalue('alternative_sessions',allChAT{36});
firingrate_analysis(nextsess{1}{1})

% SLEEP
% High firing seem to be correlated with wake-ups; FR lowest, when spectrum
% is 'empty'. So there seem to be a sleep correlation, but somewhat
% complex; FR-speed correlation also does not work out.
firingrate_analysis_sleep(nextsess{1}{2})

% No noise session
% FR does not seem to follow behavior.
firingrate_analysis(nextsess{1}{3})

% Go tone changed to 1 kHz
% All headfixed. FR maybe follows FA rate.
firingrate_analysis(nextsess{1}{4})