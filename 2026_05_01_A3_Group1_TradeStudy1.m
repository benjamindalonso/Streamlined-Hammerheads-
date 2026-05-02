clearvars
clc
close all
%% Values to consider for independent variables
AR = [1 1.5 2 2.5 3 4];
AR_continous = linspace(1,4,50);
sweepLE = [30 40 50 55 60];
sweepLE_continous = linspace(30,60,50);

% Static wing parameters for weight equation (Raymer 15.1)
Kdw = 1; % multiplier for delta wing
Kvs = 1; % multiplier for variable sweep
Wdg = 58619; % Design gross weight lb
Nz = 7; % Load factor
Sw = 600; % wing area ft^2
t_c = 0.1; % thickness ratio at root
taper = 0.24;
Scsw = 146.77 + 80.026; % Ctrl Surf Area ft^2

% Constant sweep lines
e_constSweep = spanEfficiency(AR_continous',sweepLE);
K_constSweep = (pi.*AR_continous'.*e_constSweep).^-1;
W_constSweep = wingWeight(Kdw,Kvs,Wdg,Nz,Sw,AR_continous',t_c,taper, ...
    sweepLE,Scsw)./1000; % 1000 lb

% Constant AR lines
e_constAR = spanEfficiency(AR,sweepLE_continous');
K_constAR = (pi.*AR.*e_constAR).^-1;
W_constAR = wingWeight(Kdw,Kvs,Wdg,Nz,Sw,AR,t_c,taper,sweepLE_continous' ...
    ,Scsw)./1000; % 1000 lb

% Our current design point
e_designPoint = spanEfficiency(2.027,47.3)
K_designPoint = (pi.*2.027.*e_designPoint).^-1
W_designPoint = wingWeight(Kdw,Kvs,Wdg,Nz,Sw,2.027,t_c,taper,47.3 ...
    ,Scsw)./1000 % 1000 lb
%%
plot(W_constSweep,K_constSweep,W_constAR,K_constAR,'--', ...
    W_designPoint,K_designPoint,'pr','MarkerSize',10,'MarkerFaceColor','r')
legend('Sweep = 30 deg','Sweep = 40 deg','Sweep = 50 deg','Sweep = 55 deg', ...
    'Sweep = 60 deg','AR = 1','AR = 1.5','AR = 2','AR = 2.5','AR = 3', ...
    'AR = 4','Design Point',Location='best')
xlabel('Wing Weight (x1000 lb)')
%ylim([0.55,1])
ylabel('Induced Drag Factor K')
pbaspect([2,1,1])
annotation('textarrow',[0.4,0.5]+0.2,[0.54-0.08,0.5-0.08],'String','AR ')
annotation('textarrow',[0.22,0.3],[0.7,0.8]-0.07,'String','Sweep')
grid on