function [y_fit_p,y_fit_Offset] = Stepwise_polyfit(x,y,Opt)

range = Opt.range;

y_fit_p = nan(size(x));
y_fit_Offset = nan(size(x));
for k = 1+range:length(y_fit_p)-range
    [p,S] = polyfit(x(k-range:k+range),y(k-range:k+range),1);
    y_fit_p(k) = p(1);
    y_fit_Offset(k) = p(2);
end
