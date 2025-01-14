clear all;
clc;

xs = linspace(0.1,0.9,10);
ys = zeros(1,10);

ps = [[xs;ys]]

tic;
poincarediskbaricenter(ps)'
toc

%pause

rs = randi(1000, 1, 10);
rs = rs ./ 1000;

%new_ps = [[rs;ys]]

%tic;
%poincarediskbaricenter(new_ps)
%toc
