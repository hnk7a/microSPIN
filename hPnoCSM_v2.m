% Load and displacement corretions on MI data
% Following description of Pathak et al. Determination of an effective
% zero-point and extraction of indentation stress–strain curves without the
% continuous stiffness measurement signal. 2009. Scripta Mat.

function [hs, Ps, in_P23, Res] = hPnoCSM_v2(ysec, xsec, Prange, N)

%% Inputs
% ysec is the P data
% xsec is the h data
% both for the elastic segment
% hrange is the range for hstar min to max given in a vector
% Prange same as hstar_r but for Pstar
% N is the total number of hstar and Pstars to try

%% Outputs
% hs is the displacment correction
% Ps is the load correction
% in_P23 is the index of the last imaginary value of P^(2/3)
% Res is the log10 value of the residual for hs and Ps

%% Make a arrays of possible h* and P* and evaluate

% hstar = linspace(hrange(1), hrange(2), N);
Pstar = linspace(Prange(1), Prange(2), N);
Res_ij = zeros(N,1);
in = zeros(N,1);
hstar = zeros(N,1);


for jj = 1:N
    Pnew = ysec - Pstar(jj);
    Pnew23 = Pnew.^(2/3);
    if any(imag(Pnew23)) == 1; % check for imaginary numbers
        in(jj) = find( imag(Pnew23)~=0, 1, 'last');
        % resize xsec and Pnew23
        Pnew23_i = Pnew23(in+1:end);
        xsec = xsec(in+1:end);
        s = mypolyfit(Pnew23_i, xsec, 1);
        hstar(jj) = s(2);
        xnew = xsec - s(2);
        xfit = s(1).*Pnew23_i;
        [~, ~, ~, AvgAbRes] = rsquare(xnew, xfit);
        Res_ij(jj) = log10(AvgAbRes);
        % not sure if this works
    else
        s = mypolyfit(Pnew23, xsec, 1);
        hstar(jj) = s(2);
        xnew = xsec - s(2);
        xfit = s(1).*Pnew23;
        [~, ~, ~, AvgAbRes] = rsquare(xnew, xfit);
        Res_ij(jj) = log10(AvgAbRes);
        in(jj) = 0;
    end
end


% for plotting
plot(Pstar,Res_ij,'-')
xlabel('P^*')
ylabel('log_{10}(Res)')
hold on

% find minimum Residual and corresponding P*,h*
[Res, ind] = min(Res_ij);
hs = hstar(ind);
Ps = Pstar(ind);
in_P23 = in(ind);

% plot answer
scatter(Ps, Res, 'ro', 'filled')
title(['h^*= ', num2str(hs),'um', ', ', 'P^*= ', num2str(Ps), 'N'])
hold off
        
    