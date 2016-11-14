%% PPI analysis programme
% Developed by Tavaninja
close all;
clear all;
%% Auto-Load Variables
% If Auto-Load fails, please comment the section and use the manual method

[filename1,pathname1] = uigetfile;
if  filename1 == 0;
    filename1 = 'C:\Users\metal\Desktop\master\BMAM21_3 beat.txt';
else
    filename1 = strcat(pathname1,filename1);
end
    startRow1 = 1;
    endRow1 = inf;
    VarName1 = importfile_beat(filename1, startRow1, endRow1);

[filename2,pathname2] = uigetfile;
if  filename2 == 0;
    filename2 = 'C:\Users\metal\Desktop\master\BMAM21_3 track.txt';
else
    filename2 = strcat(pathname2,filename2);
end
    startRow2 = 2;
    endRow2 = inf;
    [eTim,sObjectsFoun] = importfile_func(filename2, startRow2, endRow2);

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
% Change with the current video exact Time Size
prompt = 'Minutes? ';
mm = input(prompt)
prompt = 'Seconds? ';
ss = input(prompt)
prompt = 'MilliSeconds? ';
ms = input(prompt)
% mm=26;      % Minutes
% ss=42;      % Seconds
% ms=966;     % MiliSeconds


% Compute Total Time
TotalTime = (mm*60)+ss+ms/1000;

% % Number of Frames Back & Forth to Compute the 
% prompt = 'Seconds Back n Forth? ';
% Secs = input(prompt)
% FramesBnF = round(size(Frames,1)*Secs/TotalTime);
FramesBnF = 20;


prompt = 'Press 1 for Area, press 2 for Analytical? ';
mode = input(prompt)

%% Computations
% % Compute Total Time
% TotalTime = (mm*60)+ss+(ms/1000);

% Real & Simulated Time from 0-to-1
A = VarName1/(TotalTime);
B = Time/Time(size(Frames,1));

% Initialize a zero Vector with the size of B
C = zeros(size(B));

%% Finding The Spikes Position
for i = 1:size(A);
    for j = 1:size(B)-1;
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
if mode == 1;
    xb = zeros(1,FramesBnF);
    xf = zeros(1,FramesBnF);
    TrapezoidF = zeros(1,size(Frames,1));
    TrapezoidB = zeros(1,size(Frames,1));
end

if mode == 2;
    xb = zeros(1,FramesBnF);
    xf = zeros(1,FramesBnF);
end

% use trapezoid method for every sound
if mode == 1;
    for k = 1:size(C);
        if C(k) == -1 || C(k) == -2;
            % After the sound points for computation
            for l = 1:FramesBnF;
                xf(l) = sObjectsFoun(k+l);
            end
            % Trapezoid computation
                    TrapezoidF(k) = trapz(xf);
            % Before the sound points for computation
            for l = 1:5;
                 xb(l) = sObjectsFoun(k-l);
            end
            % Trapezoid computation
                TrapezoidB(k) = trapz(xb);
        else
            % Zero everywhere else so we dont lose the frames and easier plot
            TrapezoidB(k) = -12;
            TrapezoidF(k) = -12;
        end
    end
    
%% Analutical method
elseif mode == 2;
    for k = 1:size(C)
        if C(k) == -1 || C(k) == -2;
            for l = 1:FramesBnF;
                xf(k,l) = sObjectsFoun(k+l);
            end
            
            for l = 1:5;
                xb(k,l) = sObjectsFoun(k-l);
            end
        end 
    end
    xf(size(C),:) = 0;
    xb(size(C),:) = 0;
            % Reduce the size of x
            xf(C == 0,:) = [];
            xb(C == 0,:) = [];
end
     
%% Plots
if mode == 1;
figure(1)
Plot of the mouse movment and sound
plot(Frames,sObjectsFoun,Frames,TrapezoidB,Frames,TrapezoidF);

figure(2)
Plot with Markers
 plot(Frames,sObjectsFoun,Frames,C,Frames,TrapezoidB,Frames,TrapezoidF);
axis tight
end

%% Vectores
if mode == 1;
TrapezoidB(TrapezoidB == -12) = [];
TrapezoidB = TrapezoidB';
TrapezoidF(TrapezoidF == -12) = [];
TrapezoidF = TrapezoidF';
end

C(C == 0) = [];
