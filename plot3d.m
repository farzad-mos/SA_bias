% load('NKG2015zt.mat')
SA_mean.nkg=griddata(nkglat,nkglon,nkg2015,SA_mean.lat,SA_mean.lon');

%% 3d

close all

plot3(SA_mean.lon(SA_mean.missionid==2&SA_mean.pas==158),SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==158),movmedian(SA_mean.deltadt(SA_mean.missionid==2&SA_mean.pas==158),0.7,'SamplePoints',SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==158)),'LineWidth',4)
hold on

plot3(SA_mean.lon(SA_mean.missionid==2&SA_mean.pas==272),SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==272),movmedian(SA_mean.deltadt(SA_mean.missionid==2&SA_mean.pas==272),0.7,'SamplePoints',SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==272)),'LineWidth',4)


ylim([60.7 62.2])





figure(2)

plot3(SA_mean.lon(SA_mean.missionid==2&SA_mean.pas==158),SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==158),SA_mean.nkg(SA_mean.missionid==2&SA_mean.pas==158),'-m','LineWidth',4)
hold on
plot3(SA_mean.lon(SA_mean.missionid==2&SA_mean.pas==272),SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==272),SA_mean.nkg(SA_mean.missionid==2&SA_mean.pas==272),'-m','LineWidth',4)
ylim([60.7 62.2])
%% ascending
close all

yyaxis left

h(1)=plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==158),movmedian(SA_mean.deltadt(SA_mean.missionid==2&SA_mean.pas==158),0.7,...
    'SamplePoints',SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==158)),'-b','LineWidth',4,'DisplayName','\DeltaDT_{SA-HDM}');
hold on

plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==272),movmedian(SA_mean.deltadt(SA_mean.missionid==2&SA_mean.pas==272),...
    0.7,'SamplePoints',SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==272)),'--b','LineWidth',4)

% plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==386),movmedian(SA_mean.deltadt(SA_mean.missionid==2&SA_mean.pas==386),...
%     0.7,'SamplePoints',SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==386)),'-.b','LineWidth',4)

set(gca,'YColor','b')

xlim([60.7 64])
ylabel('\DeltaDT_{SA-HDM} [cm]','FontSize',20,'FontWeight','bold');

yyaxis right
set(gca,'YColor','m')

h(2)=plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==158),SA_mean.nkg(SA_mean.missionid==2&SA_mean.pas==158),'-m','LineWidth',4,'DisplayName','N');
hold on
plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==272),SA_mean.nkg(SA_mean.missionid==2&SA_mean.pas==272),'--m','LineWidth',4)
% plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==386),SA_mean.nkg(SA_mean.missionid==2&SA_mean.pas==386),'-.m','LineWidth',4)


xlim([60.7 64])

ylabel('Geoid Height (N) [m]','FontSize',20,'FontWeight','bold');
xlabel('Latitude [°]','FontSize',20,'FontWeight','bold');

legend(h(1:2))
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=27; ax.FontWeight='Bold'; grid on; ax.FontName='Times New Roman';
set(gca,'fontname','Times New Roman','FontSize',18);
pbaspect([1 .3 1])

%% descending

close all


yyaxis left

h(1)=plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==397&SA_mean.lat>=61.55),movmedian(SA_mean.deltadt(SA_mean.missionid==2&SA_mean.pas==397&SA_mean.lat>=61.55),0.7,...
    'SamplePoints',SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==397&SA_mean.lat>=61.55)),'-b','LineWidth',4,'DisplayName','\DeltaDT_{SA-HDM}');
hold on

plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==283),movmedian(SA_mean.deltadt(SA_mean.missionid==2&SA_mean.pas==283),...
    0.7,'SamplePoints',SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==283)),'--b','LineWidth',4)

plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==169),movmedian(SA_mean.deltadt(SA_mean.missionid==2&SA_mean.pas==169),...
    0.7,'SamplePoints',SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==169)),'-.b','LineWidth',4)

ylim([-4 8])
xlim([60.7 63.5])
ylabel('\DeltaDT_{SA-HDM} [cm]','FontSize',20,'FontWeight','bold');
set(gca,'YColor','b')

yyaxis right
set(gca,'YColor','m')

h(2)=plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==397&SA_mean.lat>=61.55),SA_mean.nkg(SA_mean.missionid==2&SA_mean.pas==397&SA_mean.lat>=61.55),'-m','LineWidth',4,'DisplayName','N');
hold on
plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==283),SA_mean.nkg(SA_mean.missionid==2&SA_mean.pas==283),'--m','LineWidth',4)
plot(SA_mean.lat(SA_mean.missionid==2&SA_mean.pas==169),SA_mean.nkg(SA_mean.missionid==2&SA_mean.pas==169),'-.m','LineWidth',4)


xlim([60.7 63.5])

ylabel('Geoid Height (N) [m]','FontSize',20,'FontWeight','bold');
xlabel('Latitude [°]','FontSize',20,'FontWeight','bold');

legend(h(1:2))
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=27; ax.FontWeight='Bold'; grid on; ax.FontName='Times New Roman';
set(gca,'fontname','Times New Roman','FontSize',18);
pbaspect([1 .3 1])
