%Smooth max
%Common x values are 0.02-3
%Needs to be robust to 265
close
x = linspace(-265, 0, 10000);
x = [x, linspace(0, 265, 10000)];

ySharp = max(x,0);
tic
ySmooth = x.^3./x.^2 .* (atan(50*x)/pi + .5);
tansoft = toc;
a = .2;
b = 1/a;
tic
% ySmooth = a*log(1+exp(b*x));
softmax = toc;
figgy = figure(46);
figure(46);
subplot(2,1,1); plot(x,ySharp); hold on; plot(x, ySmooth);

subplot(2,1,2); plot(x(1:end-1), diff(ySmooth)./diff(x));

err =ySharp-ySmooth;
errpercent = err./x * 100;
figure; plot(x, abs(errpercent));