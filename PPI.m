%% PPI analysis programme
% Developed by Tavaninja
%% Auto-Load Variables
% % If Auto-Load fails, please comment the section and use the manual method
% filename1 = 'C:\Users\metal\OneDrive\Documents\Machine_Vision\16 11 2014 BMAM21_1_BEAT.txt';
% startRow1 = 1;
% endRow1 = inf;
% VarName1 = importfile_beat(filename1, startRow1, endRow1);
% 
% filename2 = 'C:\Users\metal\OneDrive\Documents\Machine_Vision\16 11 2014 BMAM21_1.txt';
% startRow2 = 2;
% endRow2 = inf;
% [eTim,sObjectsFoun] = importfile_func(filename2, startRow2, endRow2);

%% Manual-Load
% Load the 3rd n 4th colums of the track.txt file exported from the 
% Open Vision Control and the first column of the beats txt into the
% Workspace 

% Remove the 1st ROW if NaN
% if any(isnan(eTim)) == 1;
%     eTim(1,:) = [];
%     SFram(1,:) = [];
%     sObjectsFoun(1,:) = [];
% end

%% Initialazation
Frames = linspace(1,size(eTim,1),size(eTim,1));
Frames = Frames';

% Initialize OVC Time
Time = eTim-eTim(1);

%% User Conditions
% % Change with the current video exact Time Size
% prompt = 'Minutes? ';
% mm = input(prompt)
% prompt = 'Seconds? ';
% ss = input(prompt)
% prompt = 'MilliSeconds? ';
% ms = input(prompt)
mm=29;      % Minutes
ss=59;      % Seconds
ms=670;     % MiliSeconds

% Number of Frames Back & Forth to Compute the 
% prompt = 'Seconds Back n Forth? ';
% Secs = input(prompt)
% FramesBnF = round(size(Frames,1)*Secs/TotalTime);
FramesBnF = 30;

%% Computations
% Compute Total Time
TotalTime = (mm*60)+ss+ms/1000;

% Real & Simulated Time from 0-to-1
A = VarName1/(1000000*TotalTime);
B = Time/Time(size(Frames,1));

% Initialize a zero Vector with the size of B
C = zeros(size(B));

%% Finding The Spikes Position
for i = 1:size(A);
    for j = 1:size(B);
        if (A(i) >= B(j)) && (A(i) <= B(j+1));
            C(j) = -1;
        end
    end
end
% In case of double sound remove the first Beat
for m = 1:size(C);
    for n = 1:50;
        if (C(m) == -1) && (C(m+n) == -1);
            C(m) = 0;
            C(m+n) = -2;
        end
    end
end

%% Trapezoid Method
% Compute the space using Trapezoid method 
% pre allocate vectores for speed
x = zeros(1,FramesBnF);
y = zeros(1,FramesBnF);
TrapezoidF = zeros(1,size(Frames,1));
TrapezoidB = zeros(1,size(Frames,1));
% use trapezoid method for every sound
for k = 1:size(C);
    if C(k) == -1 || C(k) == -2;
        % After the sound points for computation
        for l = 1:FramesBnF;
            y(l) = k+l;
            x(l) = B(k+l);
        end
        % Trapezoid computation
        TrapezoidF(k) = trapz(x,y);
        % Before the sound points for computation
        for l = 1:FramesBnF;
            y(l) = k-l;
            x(l) = B(k-l);
        end
        % Trapezoid computation
        TrapezoidB(k) = trapz(x,y);
    else
        % Zero everywhere else so we dont lose the frames and easier plot
        TrapezoidB(k) = 0;
        TrapezoidF(k) = 0;
    end
end

%% Plots
% Plot of the mouse movment and sound
plot(Frames,sObjectsFoun,Frames,TrapezoidB,Frames,TrapezoidF);
% Plot with Markers
% plot(Frames,sObjectsFoun,Frames,C,Frames,TrapezoidB,Frames,TrapezoidF);
axis tight

%% Vectores
TrapezoidB(TrapezoidB == 0) = [];
TrapezoidF(TrapezoidF == 0) = [];

C(C == 0) = [];