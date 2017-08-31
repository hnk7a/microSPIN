
function [data] = loadunload(fname, load_start, load_step, radius, hold_time)

fz=18; % font size for plots

% Load Data

% for an excel file
% rawdata = xlsread(fname);

% for a .tra file
rawdata = dlmread(fname,',',2,0);

data = struct( ...
    'time', rawdata(1:end,1), ...
    'depth', rawdata(1:end,2), ...
    'force', rawdata(1:end,3) ); % units s, um, N

% % if you want to limit the data imported
% limit = find(data.force > 100, 1);
% data.time(limit:end)=[];
% data.depth(limit:end)=[];
% data.force(limit:end)=[];

figure(4)
plot(data.time, data.force, 'b.'); % rawdata plot



set(gca,'FontSize',fz)
% http://tonyfast.com/nsf-goali/2014/08/29/Peak-Finding/
% Design Filter
% all 1's with a 0 in the middle

% radius ;  % filter radius, adjust for noise
            % i.e. we don't want peak-valley-peak showing up in the same peak
            % 2x peak in the same peak is okay because we can fix that
            % later

filter = ones(radius*2+1, 1);
filter(radius + 1) = 0;

% Dilate Data Using Filter and Find Peaks

data = setfield( data, 'dilate', imdilate( data.force, filter) );

% peaks exist where data > data.dilate
time_up = find( data.force >= data.dilate );

% Erode Data Using Filter and Find Valleys

data = setfield( data, 'erode', imerode( data.force, filter ) );

% valleys exist where data > data.erode
time_dn = find( data.force <= data.erode );

% Check peaks and valleys

load_max = max(data.force);
load_max_adj = load_max - (load_start - load_step);
pvs = round(load_max_adj/load_step); % expected number of peaks and valleys
% hold_time; % not seconds, more like data points in a hold or a window for only keeping peak/valley

% get rid of peaks/valleys in initial load approach
dn = find(data.force(time_dn) < 0.1);
dny = isempty(dn);
if dny == 0;
    time_dn(dn) = []; %get rid of valleys below 0.5 N 
end
up = find(data.force(time_up) < 0.5);
upy = isempty(up);
if upy == 0;
    time_up(up) = [];
end

% get rid of last peak if it doesn't have a valley after it
if time_up(end) > time_dn(end)
    time_up(end) = [];
end

% get rid of multiple peaks or valleys in the same hold
pmv = length(time_up) - length(time_dn); % should equal 0
if pmv ~= 0 || length(time_up) ~= pvs;

    % check for extra valleys
    A = squareform(pdist(time_dn));
    rvi = find(diag(A,1) <= hold_time);
    
    % remove extra valleys
    time_dn(rvi) = [];
        
    % check for extra peaks
    B = squareform(pdist(time_up));
    rpi = find(diag(B,1) <= hold_time);
    
    % remove extra peaks
    time_up(rpi) = [];
    
end

data.time_up = time_up;
data.time_dn = time_dn;
data.peaks = data.force( time_up );
data.valleys = data.force( time_dn );

data.file_name = fname;
data.filter_radius = radius;
data.filter_gap = hold_time;
data.first_unload = load_start;
data.load_step = load_step;


hold on
plot(data.time(time_up), data.peaks,'g*')
plot(data.time(time_dn), data.valleys, 'r^')

plot(data.time(20:120), data.force(20:120), 'r.');

coefficients = polyfit(data.time(20:120), data.force(20:120), 1);
% Now get the slope, which is the first coefficient in the array:
slope = coefficients(1)

pmv = length(time_up) - length(time_dn); % should equal 0
if pmv ~= 0 || length(time_up) ~= pvs;
    message = 'wrong peaks and valleys or load step'
    no_peaks = length(time_up)
    no_valleys = length(time_dn)
    expected = pvs
    return
end

title('Load [N] vs Time', 'FontSize',fz);
h=legend('all','peaks','valleys');
set(h,'Location','NorthWest')
hold off

% print the first ten peak and valley indices
num_print = 15;
if num_print > pvs;
    num_print = pvs;
end
valley_choices = time_dn(1:num_print)
peak_choices = time_up(1:num_print)

% going peak to peak uses an unload and reload
% going valley to peak uses one reload
% going peak to valley uses an unload