function Results = MicroISS_Ph_v2(data, start, stop, Ei, vi, vs, R, unload_percent, Es_m, testno)

fz = 8; % font size for plots,18

x = data.depth; % um
y = data.force; % N
time_up = data.time_up;
time_dn = data.time_dn;


%Record these parameters in Results
Results.start = start;
Results.stop = stop;
Results.Ei = Ei;
Results.vi = vi;
Results.vs = vs;
Results.R = R;
Results.unload_percent = unload_percent;
Results.Es_m = Es_m;

%% Load displacement plot (upper left) %%

figure()
subplot(2,2,1)
hold on
plot(x, y, 'b-', x(start:stop), y(start:stop), 'r*') %% red segment in load-displacement curve
plot(x(time_up), y(time_up),'g*') %% loading points
plot(x(time_dn), y(time_dn), 'r^') %% unloading points
set(gca,'FontSize',fz)
% title('Load [N] vs Displ. [um]','FontSize',fz);

title([['Load [N] vs Displ. [um]'],[[', ' num2str(start)] ['-' num2str(stop)] [', test-' num2str(testno)]] ] ,'FontSize',fz);
h=legend('all','elastic','peaks','valleys');
set(h,'Location','SouthEast')

% zoomed in
figure(2)
plot(x(1:stop+50),y(1:stop+50),'b.', x(start:stop),y(start:stop),'r*')
set(gca,'FontSize',fz)
title('Load [N] vs Displ. [um] Zoomed In','FontSize',fz);
h=legend('all','elastic');
set(h,'Location','SouthEast')

figure(1)
%% Zero Point and Modulus Determination %%

% - subscript h for all data up to the first unload 
%  ( first positive displacement value to time_up(1) )
% - subscript sec for the data of the elastic selection 
%  ( start to stop )

x_sec = x(start:stop); % zero point fit and elastic selection
y_sec = y(start:stop);

Prange = [0 4]; % Pstar range
N = 800; % number of data points for each

subplot(2,2,2)
[hstar, Pstar, in_P23, Res] = hPnoCSM_v2(y_sec, x_sec, Prange, N);

Results.hstar = hstar;
Results.Pstar = Pstar;
Results.start2 = start + in_P23;
start = start + in_P23;
y_new = y - Pstar; % apply load correction to all data
y_nsec = y_new(start:stop);
ysec23 = y_nsec.^(2/3);

x_new = x - hstar; % apply displacement correction to all data
x_nsec = x_new(start:stop);

xneg = sign(x_nsec);
if any(xneg(:)==-1) == 1;
    warning = ['corrected displacement contains negative numbers, ',...
               'pick a higher start point']
end
a_sec = (R.*x_nsec).^(0.5); % contact radius of elastic selection, um
% ae = max(a_sec); % contact radius at selected elastic limit, um
% he = max(x_nsec); % indentatin depth at elastic limit

s = mypolyfit(ysec23, x_nsec, 1);
xfit = s(1).*ysec23 + s(2);
R2_0 = rsquare(x_nsec, xfit);
hre = s(2); % residual height of elastic segment, um
k_0 = s(1); % um / N^2/3 

subplot(2,2,3)
plot(x_nsec, ysec23, '.', xfit, ysec23, '-')
set(gca,'FontSize',fz)
axis([0.9*min(x_nsec) 1.1*max(x_nsec) 0 max(ysec23)]);
title('h [um] vs p(2/3) [N]','FontSize',fz)
text(min(x_nsec), 0.10*max(ysec23), ['h_r=',num2str(hre), ', ', 'R2=',num2str(R2_0)],'FontSize',fz);

Eeff = 3/4*R^(-0.5)*k_0^(-3/2) * 10^3; % Effective Modulus, N/um^2*10^3 (GPa)
Es =(1/Eeff - (1-vi^2)/Ei)^-1*(1-vs^2); % Sample Modulus, GPa
Eind =(1/Eeff - (1-vi^2)/Ei)^-1; % Indentation Modulus, GPa

% data to first unload
in_neg = find(x_new < 0, 1, 'last');
if isempty(in_neg)==1;
    in_pos = 1;
else
    in_pos = in_neg + 1;
end

xh = x_new(in_pos:time_up(1)); % displacement during elastic loading, um
yh = y_new(in_pos:time_up(1)); % load during elastic loading, N
ah = (R.*x_new(in_pos:time_up(1))).^(0.5); % contact radius during elastic loading, um

%% Calculate the unloading stiffness and contact area
y23 = y_new.^(2/3);
ii = find(time_up >= stop, 1); % index of first unloading peak in plastic regime

unload_elastic_x = [];
unload_elastic_y = [];

for i = 1:length(time_up)- ii + 1;
       
    v = time_dn(ii+i-1);
    p = time_up(ii+i-1);
    
    
    if i == length(time_up)- ii + 1;
    else
    p1 = time_up(ii+i);
    end
    
%     a = p+50;
%     
%     b = a + 100;

    ai = find(y_new(p:v)< unload_percent(1)*y_new(p),1); % i.e. limit fit to data <95% Peak Load
    a = ai + p -1; % global index

    

    
    if unload_percent(2) > 0;
        bi = find(y_new(p:v)< unload_percent(2)*y_new(p),1); % i.e. limit fit to data >60% Peak Load
        b = bi + p -1; % global index
    end
    if unload_percent(2) == 0;
        b = v; % end of unloading even if partial unloading was used during test
    end
    
    di = find(y_new(v:p1)< unload_percent(1)*y_new(p),1,'last'); % i.e. limit fit to data <95% Peak Load
    d = di + v +1; % global index    
    
    ci = find(y_new(v:p1)< unload_percent(2)*y_new(p),1,'last'); % i.e. limit fit to data >60% Peak Load
    c = ci + v +1; % global index   

%     ai = find(y_new(p:v)< unload_percent(1)*y_new(p),1); % i.e. limit fit to data <95% Peak Load
%     a = ai + p -1; % global index
% 
% 
%     
%     if unload_percent(2) > 0;
%         bi = find(y_new(p:v)< unload_percent(2)*y_new(p),1); % i.e. limit fit to data >60% Peak Load
%         b = bi + p -1; % global index
%     end
%     if unload_percent(2) == 0;
%         b = v; % end of unloading even if partial unloading was used during test
%     end

    su = mypolyfit([y23(a:b); y23(c:d)], [x_new(a:b); x_new(c:d)], 1);
%     xfit = su(1)*y23(a:b) + su(2);
%     R2u(i,1) = rsquare(x_new(a:b), xfit);
    
    hr(i,1) = su(2); % um
    k(i,1) = su(1); % um / N^(2/3)
    hmax(i,1) = x_new(p) ;
    pmax(i,1) = y_new(p);
    Reff(i,1) = ( 4/3*k(i)^(3/2)* (Eeff*10^-3) )^-2; % um
    ap(i,1) = sqrt(Reff(i)*( hmax(i) - hr(i) )); % um  
    S(i,1) = 3/2*k(i)^(-3/2)*(hmax(i)-hr(i))^0.5; % Stiffness N/um
%     apS(i,1) = S(i)/(2*Eeff*10^-3); % this is equal to ap

% 
    unload_elastic_x = [unload_elastic_x  x_new(a:b)' x_new(c:d)'];
    unload_elastic_y = [unload_elastic_y  y_new(a:b)' y_new(c:d)'];
    
%         unload_elastic_x = [unload_elastic_x  x_new(c:d)'];
%     unload_elastic_y = [unload_elastic_y  y_new(c:d)'];
end

figure(3)
hold on

plot(x_new,y_new,'b-');
% plot(x, y, 'b-', x(start:stop), y(start:stop), 'r*') %% red segment in load-displacement curve
plot(x(time_up)-hstar, y(time_up)-Pstar,'g*') %% loading points
plot(x(time_dn)-hstar, y(time_dn)-Pstar, 'r^') %% unloading points
% set(gca,'FontSize',fz)
% % title('Load [N] vs Displ. [um]','FontSize',fz);
% 
title('(Corrected) Load vs Displacement' ,'FontSize',fz);
xlabel('Displacement (um)');
ylabel('Load (N)');
% h=legend('all','elastic','peaks','valleys');
% set(h,'Location','SouthEast')


hold on
plot(unload_elastic_x,unload_elastic_y,'k.');

h = legend('all','peaks','valleys','elastic unloading');
set(h,'Location','SouthEast')
% hold off

%% Calculate Stress and Strain

figure(1)
% the indentation depth needs to be corrected for the displacement of the
% indenter tip in order to find the sample strain and compute indentation
% properties
% subscripts: 
% h denotes data to first unload, 
% sec denotes elastic selection,
% p denotes post elastic data

x_sec_tip = 3/4*(1-vi^2).*y_nsec./(Ei.*a_sec).*10^3; % displacement of the indenter tip, um
x_sec_sample = x_nsec - x_sec_tip; % sample displacement, um

xp_tip = 3/4*(1-vi^2).*pmax./(Ei.*ap).*10^3; % displacement of the indenter tip, um
xp_sample = hmax - xp_tip; % sample displacement, um

a_secmm = a_sec.*10^-3; % contact radius in mm
apmm = ap.*10^-3; % contact radius in mm

stress_sec = y_nsec./ (pi.*a_secmm.^2); % elastic selection indentation stress, MPa
strain_sec = 4/(3*pi).* x_sec_sample./a_sec; % ealstic regime indentation strain

stress_p = pmax./(pi.*apmm.^2); % Indentation Stress, N/mm^2 (MPa)
strain_p = 4/(3*pi).*xp_sample./ap; %Indentatin Strain, um/um

Stress = [stress_sec; stress_p]; % MPa
Strain = [strain_sec; strain_p];
a = [a_sec; ap]; % um

% look at the all data up to the first unload
xh_tip = 3/4*(1-vi^2).*yh./(Ei.*ah).*10^3; % displacement of the indenter tip, um
xh_sample = xh - xh_tip; % sample displacement, um

ahmm = ah.*10^-3; % contact radius in mm
stress_h = yh./ (pi.*ahmm.^2); % elastic selection indentation stress, MPa
strain_h = 4/(3*pi).* xh_sample./ah; % ealstic regime indentation strain

%% Calculate indentation yield strength
strain_offset = 0.002;
f1 = Eind*10^3.*(Strain - strain_offset);
f2 = Stress;
YS_up = find(f2-f1 < 0, 1, 'first');

if isempty(YS_up)==1;
    Yield_Stress = 0;
    Yield_Strain = 0;
    Yield_a = 0;
else
    YS_low = YS_up -1;
    YS_y = [Stress(YS_low), Stress(YS_up)];
    YS_x = [Strain(YS_low), Strain(YS_up)];
    p2 = mypolyfit(YS_x, YS_y, 1);
    p1 = [Eind*10^3, -1*strain_offset*Eind*10^3];
    Yield_Strain = fzero(@(x) polyval(p1-p2,x),1);
    Yield_Stress = polyval(p2, Yield_Strain);
    
    % contact radius at yield
    % this is a little buggy
    Yield_a = a(YS_low) + (Yield_Stress - YS_y(1)) / (YS_y(2)-YS_y(1))/(a(YS_up)-a(YS_low)); % um
end

%% Calculate Hardening
% simple case: linear fit from Yind to last data point
H = mypolyfit(Strain(YS_up:end), Stress(YS_up:end), 1);
Hfit = H(1).*Strain(YS_up:end)+H(2);
R2H = rsquare(Stress(YS_up:end), Hfit);

Results.Hardening = H(1);
Results.HardeningFitR2 = R2H;
%% Stress-Strain Plot

Mstrain = [0 Yield_Strain];
Mstress = Eind*10^3.*Mstrain;

Ostrain = [0 Strain(end)];
Ostress = Eind*10^3.*(Ostrain - strain_offset);

if Es_m == 0
    
    subplot(2,2,4)
    hold on
    
%     plot(strain_h, stress_h, '*','MarkerEdgeColor',[0.5 0.5 0.5])
    plot(strain_sec, stress_sec, 'r.', strain_p, stress_p, 'b-*')
    set(gca,'FontSize',fz)
    
    plot(Mstrain, Mstress, 'k-', 'LineWidth',1)
    plot(Ostrain, Ostress, 'k--', 'LineWidth',2)
    plot(Yield_Strain, Yield_Stress,'g.','MarkerSize',25)
    
    plot(Strain(YS_up:end), Hfit,'r--', 'LineWidth',1);
    
    title(['Stress [MPa] vs Strain: ',['Es=',num2str(Es)], ', ', ['YS=', num2str(Yield_Stress)] , ', ', ['a=', num2str(Yield_a)]] ,'FontSize',fz)
    % use this when you take out the contact radius at yield calculation
    % title(['Stress [MPa] vs Strain: ',['Es=',num2str(Es)], ', ', ['YS=', num2str(Yield_Stress)] ], 'FontSize',fz)
    
    h=legend('data to 1^s^t unload', 'elastic', 'post-elastic', 'modulus', 'offset', 'yield');
    % h=legend('elastic', 'post-elastic', 'modulus', 'offset', 'yield');
    set(h,'Location', 'SouthEast')
    axis([0 1.1*max(Strain), 0 1.1*max(Stress)]);
    hold off
    
    % Results into Structure
    Results.Stress = Stress;
    Results.Strain = Strain;
    Results.a = a;
    
    Results.stress_e = stress_sec;
    Results.strain_e = strain_sec;
    Results.Esample = Es;
    Results.Yield_Strength = Yield_Stress;
    Results.Yield_Strain = Yield_Strain;
    Results.Yield_a = Yield_a;
end

%% Calculation and Plot for an Assumed Modulus %%
% Es_m assumed sample modulus, GPa
% look at things if you assume a modulus and apply the zero point
% correction, subscript m denotes assumed Moudlus was used

if Es_m > 0
    
    Eeff_m = ((1-vi^2)/Ei + (1-vs^2)/Es_m)^-1; % Effective Modulus, N/um^2*10^3 (GPa)
    Eind_m =(1/Eeff_m - (1-vi^2)/Ei)^-1; % Indentation Modulus, GPa
    
    xx = x_new; % uses displacement correction
    % xx = x; % no displ. correction
    
    % the displacement correction has no effect on a k(i) or hmax(i)-hr(i),
    % however we need to include any unloads we've skipped
    
    for i = 1:length(time_up);
        
        vm = time_dn(i);
        pm = time_up(i);
        
        aim = find(y_new(pm:vm)< unload_percent(1)*y_new(pm),1); % i.e. limit fit to data <95% Peak Load
        am = aim + pm -1; % global index
        if unload_percent(2) > 0;
            bim = find(y_new(pm:vm)< unload_percent(2)*y_new(pm),1); % i.e. limit fit to data >60% Peak Load
            bm = bim + pm -1; % global index
        end
        if unload_percent(2) == 0;
            bm = vm; % end of unloading even if partial unloading was used during test
        end
        
        su_m = mypolyfit(y23(am:bm), xx(am:bm), 1);
        xfit = su_m(1)*y23(am:bm) + sum(2);
        R2u_m(i,1) = rsquare(xx(am:bm), xfit);
        
        hr_m(i,1) = su_m(2); % um
        k_m(i,1) = su_m(1); % um / N^(2/3)
        hmax_m(i,1) = xx(pm);
        pmax_m(i,1) = y_new(pm);
        Reff_m(i,1) = ( 4/3*k_m(i)^(3/2)* (Eeff_m*10^-3) )^-2; % um
        a_m(i,1) = sqrt(Reff_m(i)*( hmax_m(i) - hr_m(i) )); % um
        Sm(i,1) = 3/2*k_m(i)^(-3/2)*(hmax_m(i)-hr_m(i))^0.5; % Stiffness N/um
        %     apS(i,1) = S(i)/(2*Eeff*10^-3); % this is equal to ap
    end
    
    xm_tip = 3/4*(1-vi^2).*pmax_m./(Ei.*a_m).*10^3; % displacement of the indenter tip, um
    xm_sample = hmax_m - xm_tip; % sample displacement, um
    
    am_mm = a_m.*10^-3; % contact radius in mm
    stress_m = pmax_m./ (pi.*am_mm.^2); % indentation stress, MPa
    strain_m = 4/(3*pi).* xm_sample./a_m; % indentation strain
    
    % figure()
    subplot(2,2,4)
    hold on
    plot(strain_h, stress_h, '*','MarkerEdgeColor',[0.5 0.5 0.5])
    plot(strain_m, stress_m, 'm*', strain_sec, stress_sec, 'r.', strain_p, stress_p, 'b-*')
    set(gca,'FontSize',fz)
    Mstress_m = Eind_m*10^3.*Mstrain;
    plot(Mstrain, Mstress_m, 'k-', 'LineWidth',1)
    title('Stress [MPa] vs Strain','FontSize',fz)
    h=legend('data to 1^s^t unload w/ 0pt','Assumed Modulus', 'elastic w/ 0pt', 'post-elastic', 'Modulus');
    set(h,'Location', 'SouthEast')
    text(0.1*max(strain_sec), 0.10*max(stress_sec), ['Es=',num2str(Es_m)],'FontSize',fz);
    axis([0 1.1*max(strain_m), 0 1.1*max(stress_m)]);
    
    % find the YS of the assumed modulus stress-strain curve
    f1 = Eind_m*10^3.*(strain_m - strain_offset);
    f2 = stress_m;
    YS_up = find(f2-f1 < 0, 1, 'first');
    if isempty(YS_up)==0 && YS_up > 1 ;
        YS_low = YS_up -1;
        YS_y = [stress_m(YS_low), stress_m(YS_up)];
        YS_x = [strain_m(YS_low), strain_m(YS_up)];
        p2 = mypolyfit(YS_x, YS_y, 1);
        p1 = [Eind_m*10^3, -1*strain_offset*Eind_m*10^3];
        Yield_Strain_m = fzero(@(x) polyval(p1-p2,x),1);
        Yield_Stress_m = polyval(p2, Yield_Strain_m);
        
        % contact radius at yield
        % this is buggy
        Yield_a_m = a_m(YS_low) + (Yield_Stress_m - YS_y(1)) / (YS_y(2)-YS_y(1))/(a_m(YS_up)-a_m(YS_low)); % um
        
    else
        % if there is no data to the left of the strian offset line, take the
        % first stress-strain point
        Yield_Strain_m = strain_m(1);
        Yield_Stress_m = stress_m(1);
        Yield_a_m = a_m(1);
    end
    
    Ostrain = [0 Yield_Strain_m];
    Ostress = Eind_m*10^3.*(Ostrain - strain_offset);
    
    plot(Ostrain, Ostress, 'k--', 'LineWidth',2)
    plot(Yield_Strain_m, Yield_Stress_m,'g.','MarkerSize',25)
    
    title(['Stress [MPa] vs Strain: ',['Es=',num2str(Es_m)], ', ', ['YS=', num2str(Yield_Stress_m)] , ', ', ['a=', num2str(Yield_a_m)]] ,'FontSize',fz)
    hold off
    
    Results.Stress = stress_m;
    Results.Strain = strain_m;
    Results.a = a_m;
    Results.Esample = Es_m;
    Results.Yield_Strength = Yield_Stress_m;
    Results.Yield_Strain = Yield_Strain_m;
    Results.Yield_a = Yield_a_m;
end
