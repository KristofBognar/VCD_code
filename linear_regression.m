function [m, Y_int, R2, Y_fit, rms] = linear_regression(X,Y, alpha, fig_flag)
% This function performs and plots a linear regression on the given dataset
% Written Dec 2010 by Cristen Adams
%
% INPUT
% X = X-vector
% Y = Y-vector
% alpha = confidence interval is equal to: 100*(1-alpha)%
% fig_flag = 0-no figure, 1-make figure
% 
% OUTPUT:  Just a figure with labels for R^2, m, and Y
%  m = [slope, error_in_slope]
%  Y_int = [y-intercep, error-in-y-int]
%  R2 = r-squared value for fit
%  Y_fit = Y values calculated by the fit

try 
    [b,bint,r,rint,stats] = regress(Y, [ones(length(X),1) X], alpha);
catch
    disp('Could not run regression.')
    m = [NaN NaN];
    Y_int = [NaN NaN];
    R2 = NaN;
    Y_fit = NaN;
    rms = NaN;
    return
end
% get the slope and y-intercept values with their ranges within confidence
% levels

% get y-values of the fit
Y_fit = polyval(flipud(b), X);

% deal with residuals
ind = find(~isnan(r));
rms = norm(r(ind)) / sqrt(length((ind)));
r_avg = mean(r(ind));
sigma = std(r(ind));
% R^2
R2 = stats(1);

% now get +/- errors
Y_err1 = abs(bint(1,1) - b(1));
Y_err2 = abs(bint(1,2) - b(1));
% if Y_err1 ~= Y_err2,
%     disp('Warning: error in Y-intercept is not symmetric.')
%     disp('I have no idea why this would happen...')
% end
Y_err = mean([Y_err1, Y_err2]);

m_err1 = abs(bint(2,1) - b(2));
m_err2 = abs(bint(2,2) - b(2));
% if m_err1 ~= m_err2,
%     disp('Warning: error in slope is not symmetric.')
%     disp('I have no idea why this would happen...')
% end
m_err = mean([m_err1, m_err2]);

%m = [bint(2,1), m_err]; % this is original code from Cristen. should not
%use bint(2,1), which is just the lower 95% bound of the slop
%Y_int = [bint(1,1), Y_err]; % this is also from Cristen's original code
m = [b(2), m_err];
Y_int = [b(1), Y_err];

X_axis = 1:length(X);
if fig_flag == 1,
    figure(2)
    subplot(2,1,1)
    hold off
    plot(X, Y, 'bo', 'markerfacecolor','m')
    hold on
    plot(X, Y_fit, 'k--')
    ylabel('Y')
    xlabel('X')
    m_str = [num2str(m(1),4) ' {\pm ' num2str(m(2),4) '}'];
    Y_str = [num2str(Y_int(1),4) ' {\pm ' num2str(Y_int(2),4) '}'];
    title({['R^2 = ' num2str(R2,4) ],...
        ['m = ' m_str], ['Y_{int} = ' Y_str]});
    subplot(2,1,2)
    hold off
    plot(X_axis, r, '*')
    hold on

    plot(X_axis, r_avg * ones(1,length(X)), 'm-', 'linewidth', 1)
    plot(X_axis, zeros(1,length(X)), 'k--', 'linewidth', 2)
    legend('Residuals', 'Avg residuals', 'Zero-line')
    %xlim([min(X), max(X)])
    xlabel('Index')
    ylabel('Y minus Y-fit')
    title(['RMS: ' num2str(rms, 4) '{\pm}' num2str(sigma,4)])
end
end
