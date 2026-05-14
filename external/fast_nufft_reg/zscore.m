function [x,x_mean,x_std]=zscore(x)

x_mean=mean(x);
x_std=std(x);
x=(x-x_mean)./x_std;

end