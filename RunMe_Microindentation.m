    %% Microindentation Stress-Strain Analysis %%
% Jordan Weaver, 7/17/2015, Georgia Tech MINED
close all
clear all
clc
 
%Load Data and Find Peaks and Valleys %%

% floc = '\\Client\E$\Jonny\10-10-2016\';
% fname = 'marsteel_12test-';
tnum = '10'; 
testno = str2num(tnum);
% fname = [floc,fname,tnum,'.xlsx']; % for excel
% fname = [floc,fname,tnum, '.tra']; % for .tra      %must also change file type in loadunload.m

fname = uigetfile('*.tra','select TRA');


load_start = 2; % A test parameter, load to for first unload
load_step = 2; % A test parameter, load increase per cycle
radius = 55; % filter radius for peak/valley finding
gap = 100; %min. data points between two peaks/valleys

Data = loadunload(fname, load_start, load_step, radius, gap);

%% Manually Remove an Unloading Cycle if needed
% ii = 1;
% data.time_up(ii) = [];
% data.time_dn(ii) = [];
% data.peaks(ii) = [];
% data.valleys(ii) = [];
%% Calculate Indentation Stress and Strain %%
close all
clc;
start = 20; % start of elastic segment selection, try to use a peak/valley
stop = 120; % end of elastic segment selection, try to use a peak/valley
Ei = 640;   % Indenter Modulus, GPa, Tungsten Carbide
vi = 0.211; % Indenter Poisson Ratio
vs = 0.3;   % Sample Poisson Ratio
R = 6350; % Indenter Radius, um, 6350
% 6350;  % "0.5 inch diameter"
% 2500; % data from Zwick Germany was w/ "5 mm diameter ball"
% 500; % "1mm diameter ball"
unload_percent = [.85 .55]; % portion of unloading curve for stiffness calc.
Es_expected = 00; % expected sample modulus in GPa for diagnostic plot, %105
Results = MicroISS_Ph_v2(Data, start, stop, Ei, vi, vs, R, unload_percent, Es_expected, testno);

% BELOW COMMENTED OUT TEXT WORKS FOR MATLAB 2014+
scrsz = get(groot,'ScreenSize'); %works with matlab 2014+
fsz = [scrsz(1)+100 scrsz(2)+100 scrsz(3)-200 scrsz(4)-200];
set(gcf,'Position', fsz)

% Put everything in a structure for saving
Analysis.Results = Results;
Analysis.Data = Data;

%% save analysis and the plot in your current directory using tnum
save(['Analysis ' tnum],'Analysis')
set(gcf,'PaperPositionMode','auto')
saveas(gcf,['ISS ' tnum], 'png')



