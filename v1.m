load('sa.mat')

% make sa data cell based on said
for q=1:3
    SA{q}=sa(sa.said==q,:);
end

clearvars sa q

%% 
% delta points 1 nuatical mile
deltalat=km2deg(nm2km(1));
deltalon=km2deg(nm2km(1))*2;

Lat1=(53.7:km2deg(nm2km(1)):65.8910)';
Lon1=(10.01:km2deg(nm2km(1))*2:30.1964);


Lon=repmat(Lon1,size(Lat1,1),1);
Lat=repmat(Lat1,1,size(Lon1,2));
% 
%% mean per pass

 for k=1:3
    meanSA{k}=table();
    H=1;
    pas=unique(SA{k}.Pass);
    for i=1:length(pas)
        temp=SA{k}(SA{k}.Pass==pas(i),:);
        
        cyc=unique(temp.Cycle);
        lat=[min(temp.Lat(temp.Cycle==cyc(1))):deltalat:...
            min(temp.Lat(temp.Cycle==cyc(1)))+(ceil((max(temp.Lat(temp.Cycle==cyc(1)))-min(temp.Lat(temp.Cycle==cyc(1))))/deltalat)*deltalat)]';
        
        for j=1:length(lat)
            lon(j,:)=mean(temp.Lon(temp.Lat>=lat(j)-deltalat/2&temp.Lat<lat(j)+deltalat/2));
        end
        
        for m=1:length(cyc)
            temp1=temp(temp.Cycle==cyc(m),:);
            
            for n=1:length(lat)
               samean=mean(temp1.DTsa(temp1.Lat>=lat(n)-(deltalat/2)&temp1.Lat<lat(n)+(deltalat/2)),'omitnan');
               hdmmean=mean(temp1.DThdm(temp1.Lat>=lat(n)-(deltalat/2)&temp1.Lat<lat(n)+(deltalat/2)),'omitnan');
                
                
                meanSA{k}.sa(H,1)=samean;
                meanSA{k}.hdm(H,1)=hdmmean;
                meanSA{k}.lat(H,1)=lat(n);
                meanSA{k}.lon(H,1)=lon(n);
                meanSA{k}.pas(H,1)=temp1.Pass(1);
                meanSA{k}.cycle(H,1)=cyc(m);
                meanSA{k}.time(H,1)=temp1.Time(1);
                H=H+1;
                clearvars hdmmean samean
            end
            clear temp1 n
        end
        clearvars lat lon temp cyc m j
    end
    clearvars i pas
end

clearvars i j k m n pas temp H temp1 lat lon samean hdmmean cyc

%% load data
load('meanSA.mat')

% remove nan values
for k=1:3
    meanSA{k} = rmmissing(meanSA{k},'DataVariables',{'sa'}); %remove NaN data
end

for k=1:3
    meanSA{k} = rmmissing(meanSA{k},'DataVariables',{'hdm'}); %remove NaN data
end
str=['JA3';'S3A';'S3B'];
%% decompose per pass

% decompositions
for k=1:3
    pas=unique(meanSA{k}.pas);
    for i=1:length(pas)
        temp=meanSA{k}(meanSA{k}.pas==pas(i),:);
        lat=unique(temp.lat);
        for j=1:length(lat)
            
            dt_sa=temp.sa(temp.lat==lat(j));
            dt_hdm=temp.hdm(temp.lat==lat(j));
            ti=temp.time(temp.lat==lat(j));
            delta_dt=dt_hdm-dt_sa;
            rmse=rms(delta_dt,'omitnan');

            
            % EMD decompostion
            [imf,residual,~] = emd(dt_sa,'MaxNumIMF',4);
            [imf2,residual2,~] = emd(dt_hdm,'MaxNumIMF',4);
            
            d_res=residual-residual2;
%           p1 = polyfit(decyear(ti),d_res,1);
%           tre = polyval(p1,decyear(ti));
            tr1=fitlm(decyear(ti),delta_dt);
            a1=tr1.Coefficients.Estimate(2);
            
%             tr2=fitlm(decyear(ti),delta_dt);
%             a2=tr1.Coefficients.Estimate(2);
            
            %EMD Figure
            figure(1)
            for m=1:size(imf,2)+2
                subplot(size(imf,2)+2,1,m)
                
                if m==1
                    
                     z1=dt_sa;y1=dt_hdm;
                     b = [ones(size(z1,1),1) z1]\y1;
                     RegressionLine = [ones(size(z1,1),1) z1]*b;
                     SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
                     SS_Y = sum((y1-mean(y1)).^2);
                     SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
                     R_squared = SS_XY/sqrt(SS_X*SS_Y);
                    
                    h(1)=plot(ti,dt_sa,'--*','color',[0.4940 0.1840 0.5560],'LineWidth',1.5,'DisplayName',...
                        strcat(str(k,:),', RMSE_{SA-HDM}= ',num2str(rmse,3),' [cm], ',' Mean_{SA-HDM}= ',num2str(mean(delta_dt),3),'[cm]'));
                    hold on
                    h(2)=plot(ti,dt_hdm,'--k','LineWidth',1.5,'DisplayName',strcat('HDM',', R_{SA-HDM}^{2}= ','  ',num2str(R_squared,3)));
                    hold off
                    xticks(datetime('01-Dec-2016'):calmonths(3):datetime('01-Jun-2019'))
                    xlim([datetime('01-Dec-2016'),datetime('01-Jun-2019')])
                    legend(h(1:2))
                    title('EMD Decomposition')

                elseif m==size(imf,2)+2
                    h(1)=plot(ti,residual,'color',[0.4940 0.1840 0.5560],'LineWidth',1.5,'DisplayName',strcat('residual tilt=',num2str(a1,2),' [mm/yr]'));
                    hold on
                    plot(ti,residual2,'-k','LineWidth',1.5)
                    xticks(datetime('01-Dec-2016'):calmonths(3):datetime('01-Jun-2019'))
                    xlim([datetime('01-Dec-2016'),datetime('01-Jun-2019')])
                    legend(h(1))
                else
                    plot(ti,imf(:,m-1),'color',[0.4940 0.1840 0.5560],'LineWidth',1.5)
                    hold on
                    plot(ti,imf2(:,m-1),'-k','LineWidth',1.5)
                    xticks(datetime('01-Dec-2016'):calmonths(3):datetime('01-Jun-2019'))                  
                    xlim([datetime('01-Dec-2016'),datetime('01-Jun-2019')])
                    
                end
                
                if m==1
                    ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
                elseif m==size(imf,2)+2
                    ylabel('Aprox.','FontSize',18,'FontWeight','bold');
                else
                    ylabel(strcat('IMF',num2str(m)-1),'FontSize',20,'FontWeight','bold');
                end
                
                if m~=size(imf,2)+2
                    xticklabels([])
                end
                ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
                xlim([min(ti)-2,max(ti)+2])
                
            end
            
            
            % EWT decomposition (Empirical wavelet transform)
            [mra,cfs]  = ewt(dt_sa,'MaxNumPeaks',4);
            [mra2,cfs2]  = ewt(dt_hdm,'MaxNumPeaks',4);
            
            bias=mean(mra2(:,4)-mra(:,4),'omitnan');

            %EWT Figure
            figure(2)
            for m=1:size(mra,2)+1
                subplot(size(mra,2)+1,1,m)
                if m==1
                    
                     z1=dt_sa;y1=dt_hdm;
                     b = [ones(size(z1,1),1) z1]\y1;
                     RegressionLine = [ones(size(z1,1),1) z1]*b;
                     SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
                     SS_Y = sum((y1-mean(y1)).^2);
                     SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
                     R_squared = SS_XY/sqrt(SS_X*SS_Y);
                    
                    plot(ti,dt_sa,'--*','color',[0.4940 0.1840 0.5560],'LineWidth',1.5,'DisplayName',...
                        strcat(str(k,:),', RMSE_{SA-HDM}= ',num2str(rmse,3),' [cm], ',' Mean_{SA-HDM}= ',num2str(mean(delta_dt),3),'[cm]'))
                    hold on
                    plot(ti,dt_hdm,'--k','LineWidth',1.5,'DisplayName',strcat('HDM',', R_{SA-HDM}^{2}= ','  ',num2str(R_squared,3)))
                    xticks(datetime('01-Dec-2016'):calmonths(3):datetime('01-Jun-2019'))
                    hold off
                    legend show
                    xlim([datetime('01-Dec-2016'),datetime('01-Jun-2019')])
                    title('EWT Decomposition')             
                elseif m==5
                     z1=mra(:,m-1);y1=mra2(:,m-1);
                     b = [ones(size(z1,1),1) z1]\y1;
                     RegressionLine = [ones(size(z1,1),1) z1]*b;
                     SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
                     SS_Y = sum((y1-mean(y1)).^2);
                     SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
                     R_squared = SS_XY/sqrt(SS_X*SS_Y);
                     
                     
                    h(1)=plot(ti,mra(:,m-1),'-','color',[0.4940 0.1840 0.5560],'LineWidth',1.5,'DisplayName',strcat('Bias_{SA-HDM}= ',num2str(bias,3),' [cm]'));
                    hold on
                    h(2)=plot(ti,mra2(:,m-1),'-k','LineWidth',1.5,'DisplayName',strcat('R^{2}= ','  ',num2str(R_squared,3)));
                    legend(h(1:2))
                    xticks(datetime('01-Dec-2016'):calmonths(3):datetime('01-Jun-2019'))
                else
                    
                     z1=mra(:,m-1);y1=mra2(:,m-1);
                     b = [ones(size(z1,1),1) z1]\y1;
                     RegressionLine = [ones(size(z1,1),1) z1]*b;
                     SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
                     SS_Y = sum((y1-mean(y1)).^2);
                     SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
                     R_squared = SS_XY/sqrt(SS_X*SS_Y);
                 
                    plot(ti,mra(:,m-1),'-','color',[0.4940 0.1840 0.5560],'LineWidth',1.5)
                    hold on
                    h(1)=plot(ti,mra2(:,m-1),'-k','LineWidth',1.5,'DisplayName',strcat('R^{2}= ','  ',num2str(R_squared,3)));
                    legend(h(1));
                    xticks(datetime('01-Dec-2016'):calmonths(3):datetime('01-Jun-2019'))
                end
                
                if m==1
                    ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
                else
                    ylabel(strcat('MRA',num2str(m)-1),'FontSize',18,'FontWeight','bold');
                end
                ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
                xlim([datetime('01-Dec-2016'),datetime('01-Jun-2019')])
                
                if m~=size(mra,2)+1
                    xticklabels([])
                end
            end
            
            % VMD decomposition
            [IMF,res] = vmd(dt_sa,'NumIMFs',4);
            [IMF2,res2] = vmd(dt_hdm,'NumIMFs',4);
            
            bias=mean(IMF2(:,4)-IMF(:,4),'omitnan');

            
            %VMD Figure
            figure(3)
            for m=1:size(IMF,2)+2
                subplot(size(IMF,2)+2,1,m)
                
                if m==1
                    
                     z1=dt_sa;y1=dt_hdm;
                     b = [ones(size(z1,1),1) z1]\y1;
                     RegressionLine = [ones(size(z1,1),1) z1]*b;
                     SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
                     SS_Y = sum((y1-mean(y1)).^2);
                     SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
                     R_squared = SS_XY/sqrt(SS_X*SS_Y);
                     
                     
                    h(1)=plot(ti,dt_sa,'--*','color',[0.4940 0.1840 0.5560],'LineWidth',1.5,'DisplayName',...
                        strcat(str(k,:),', RMSE_{SA-HDM}= ',num2str(rmse,3),' [cm], ',' Mean_{SA-HDM}= ',num2str(mean(delta_dt),3),'[cm]'));
                    hold on
                    h(2)=plot(ti,dt_hdm,'--k','LineWidth',1.5,'DisplayName',strcat('HDM',', R_{SA-HDM}^{2}= ','  ',num2str(R_squared,3)));
                    hold off
                    title('VMD Decomposition')             
                    legend(h(1:2)) 
                    xticks(datetime('01-Dec-2016'):calmonths(3):datetime('01-Jun-2019'))
                    xlim([datetime('01-Dec-2016'),datetime('01-Jun-2019')])
                    
                elseif m==size(IMF,2)+2
                    
                    plot(ti,res,'color',[0.4940 0.1840 0.5560],'LineWidth',1.5)
                    hold on
                    plot(ti,res2,'-k','LineWidth',1.5)
                    xlim([datetime('01-Dec-2016'),datetime('01-Jun-2019')])
                    xticks(datetime('01-Dec-2016'):calmonths(3):datetime('01-Jun-2019'))
                    
                elseif m==size(IMF,2)+1
                     z1=IMF(:,m-1);y1=IMF2(:,m-1);
                     b = [ones(size(z1,1),1) z1]\y1;
                     RegressionLine = [ones(size(z1,1),1) z1]*b;
                     SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
                     SS_Y = sum((y1-mean(y1)).^2);
                     SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
                     R_squared = SS_XY/sqrt(SS_X*SS_Y);
                     
                     
                    h(1)=plot(ti,IMF(:,m-1),'color',[0.4940 0.1840 0.5560],'LineWidth',1.5,'DisplayName',strcat('Bias_{SA-HDM}= ',num2str(bias,3),' [cm]'));
                    hold on
                    h(2)=plot(ti,IMF2(:,m-1),'-k','LineWidth',1.5,'DisplayName',strcat('R^{2}= ','  ',num2str(R_squared,3)));
                    xlim([datetime('01-Dec-2016'),datetime('01-Jun-2019')])
                    legend(h(1:2))
                    xticks(datetime('01-Dec-2016'):calmonths(3):datetime('01-Jun-2019'))
                    
                else
                    z1=IMF(:,m-1);y1=IMF2(:,m-1);
                    b = [ones(size(z1,1),1) z1]\y1;
                    RegressionLine = [ones(size(z1,1),1) z1]*b;
                    SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
                    SS_Y = sum((y1-mean(y1)).^2);
                    SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
                    R_squared = SS_XY/sqrt(SS_X*SS_Y);
                    
                    
                    plot(ti,IMF(:,m-1),'color',[0.4940 0.1840 0.5560],'LineWidth',1.5)
                    hold on
                    h(1)=plot(ti,IMF2(:,m-1),'-k','LineWidth',1.5,'DisplayName',strcat('R^{2}= ','  ',num2str(R_squared,3)));
                    legend(h(1))
                    xlim([datetime('01-Dec-2016'),datetime('01-Jun-2019')])
                    xticks(datetime('01-Dec-2016'):calmonths(3):datetime('01-Jun-2019'))
                    
                end
                
                if m==1
                    ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
                elseif m==size(IMF,2)+2
                    ylabel('Aprox.','FontSize',18,'FontWeight','bold');
                else
                    ylabel(strcat('IMF',num2str(m)-1),'FontSize',20,'FontWeight','bold');
                end
                
                if m~=size(IMF,2)+2
                    xticklabels([])
                end
                
                
                ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
                xlim([min(ti)-2,max(ti)+2])
                
            end
            
            % Wavelet decomposition
            [c,l] = wavedec(dt_sa,6,'sym8');
            approx = appcoef(c,l,'sym8');
            cd = detcoef(c,l,[1 2 3 4 5 6]);
            
            [c2,l2] = wavedec(dt_hdm,6,'sym8');
            approx1 = appcoef(c2,l2,'sym8');
            cd1 = detcoef(c2,l2,[1 2 3 4 5 6]);
            
            figure(4)
            for m=1:length(cd)+2
                subplot(length(cd)+2,1,m)
                if m==1
                    plot(ti,dt_sa,'--*','color',[0.4940 0.1840 0.5560],'LineWidth',1)
                    hold on
                    plot(ti,dt_hdm,'--k')
                    ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
                    title('Wavelet Decomposition')             
                    hold off
                    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
                    xlim([min(ti)-2,max(ti)+2])
                    
                elseif m==length(cd)+2
                    plot(approx,'-','color',[0.4940 0.1840 0.5560],'LineWidth',1)
                    hold on
                    plot(approx1,'-k')
                    ylabel('Approx. ','FontSize',18,'FontWeight','bold');
                    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
                else
                    plot(cd{m-1},'-','color',[0.4940 0.1840 0.5560],'LineWidth',1)
                    hold on
                    plot(cd1{m-1},'-k')
                    ylabel(strcat('Detail ',num2str(m)-1),'FontSize',18,'FontWeight','bold');
                    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
                end
            end
            
            
            
            %         b = [ones(size(z1,1),1) z1]\y1;
            %         RegressionLine = [ones(size(z1,1),1) z1]*b;
            %         SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
            %         SS_Y = sum((y1-mean(y1)).^2);
            %         SS_XY = sum((RegressionLine-mean(RegressionLine)).*(y1-mean(y1)));
            %         R_squared = SS_XY/sqrt(SS_X*SS_Y);
            %
            
            
        end
    end
end

%% plot one point location

for n=1:length(lat)
lon(n,1)=unique(temp.lon(temp.lat==lat(n)));
end

latlim = [53 67];
lonlim = [10 31];
ax = usamap(latlim, lonlim);
hold on
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
setm(ax, 'FFaceColor', [1 1 1])
setm(gca,'FLineWidth',5,'Grid','on','FontSize',18,'fontweight','bold')
plotm(lat,lon,'-r','Linewidth',2);
plotm(lat(j),lon(j),'ob','MarkerSize',8);
plotm(lat(j),lon(j),'*b','MarkerSize',6);

%% tilt and bias                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
for k=1:3
    H=1;
    decom{k}=table();
    pas=unique(meanSA{k}.pas);
    for i=1:length(pas)
        temp=meanSA{k}(meanSA{k}.pas==pas(i),:);
        lat=unique(temp.lat);
        for j=1:length(lat)
            
            dt_sa=temp.sa(temp.lat==lat(j));
            if length(dt_sa)>5
            dt_hdm=temp.hdm(temp.lat==lat(j));
            ti=temp.time(temp.lat==lat(j));
            delta_dt=dt_hdm-dt_sa;
            decom{k}.lon(H,1)=unique(temp.lon(temp.lat==lat(j)));
            decom{k}.lat(H,1)=lat(j);
%             rmse=rms(delta_dt,'omitnan');

            
            % EMD decompostion
            [imf,residual,~] = emd(dt_sa,'MaxNumIMF',4);
            [imf2,residual2,~] = emd(dt_hdm,'MaxNumIMF',4);
            
            d_res=residual-residual2;
            tr1=fitlm(decyear(ti),d_res);
            decom{k}.tilt(H,1)=tr1.Coefficients.Estimate(2);
            
            % VMD decomposition
            [IMF,res] = vmd(dt_sa,'NumIMFs',4);
            [IMF2,res2] = vmd(dt_hdm,'NumIMFs',4);
            
            decom{k}.bias(H,1)=mean(IMF2(:,4)-IMF(:,4),'omitnan'); 
            H=H+1;
            end
            clearvars dt_sa dt_hdm ti delta_dt imf imf2 d_res tr1 IMF IMF2 residual residual2 res res2  
        end
        clearvars temp lat
    end
end

%% plot decom
for k=1:3
figure,
latlim = [53 67];
lonlim = [10 31];
ax = usamap(latlim, lonlim);
hold on
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
setm(ax, 'FFaceColor', [1 1 1])
setm(gca,'FLineWidth',5,'Grid','on','FontSize',18,'fontweight','bold')
scatterm(decom{k}.lat,decom{k}.lon,30,decom{k}.bias,'filled','o');
c=colorbar; 
c.Label.String = 'Bias_S_A_-_H_D_M (VMD) [cm]';
caxis([0 20]); %fix the bar: SA

set(c,'position',[.71 .2 .015 .7],'FontSize',18,'FontWeight','bold')


title(strcat(str(k,:),': ','Bias_S_A_-_H_D_M','Mean=',num2str(mean(decom{k}.bias),2),' [cm]'),'FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; ax.FontName='Times New Roman';
xLoc =1.3001e+05;
yLoc =7.0876e+06;
scaleruler('Units', 'km', 'RulerStyle', 'patches', ...
'XLoc', xLoc, 'YLoc', yLoc,"FontSize",10);

end


for k=1:3
figure,
latlim = [53 67];
lonlim = [10 31];
ax = usamap(latlim, lonlim);
hold on
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
setm(ax, 'FFaceColor', [1 1 1])
setm(gca,'FLineWidth',5,'Grid','on','FontSize',18,'fontweight','bold')
scatterm(decom{k}.lat,decom{k}.lon,30,decom{k}.tilt,'filled','o');
c=colorbar; 
c.Label.String = 'tilt (EMD) [cm]';
caxis([-10 10]); %fix the bar: SA


% colormap 
newmap=colormap(bluewhitered(60)); 
colbarNo = linspace(-40,40,size(newmap,1))';
newmap(abs(colbarNo) < 2,:) = repmat([1 1 1],sum(abs(colbarNo) < 2),1);
colormap(ax,newmap)




set(c,'position',[.71 .2 .015 .7],'FontSize',18,'FontWeight','bold')


title(strcat(str(k,:),': ','tilt','Mean=',num2str(mean(decom{k}.tilt),2)),'FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; ax.FontName='Times New Roman';
xLoc =1.3001e+05;
yLoc =7.0876e+06;
scaleruler('Units', 'km', 'RulerStyle', 'patches', ...
'XLoc', xLoc, 'YLoc', yLoc,"FontSize",10);

end

%% bias by mean
sabias=table();
H=1;
for k=1:3
    pas=unique(meanSA{k}.pas);
    for i=1:length(pas)
        temp=meanSA{k}(meanSA{k}.pas==pas(i),:);
        lat=unique(temp.lat);
        lon=unique(temp.lon);

        for j=1:length(lat)
            sabias.bias(H)=mean(temp.hdm(temp.lat==lat(j))-temp.sa(temp.lat==lat(j)),'omitnan');
            sabias.lat(H)=lat(j);
            sabias.lon(H)=unique(temp.lon(temp.lat==lat(j)));
            sabias.pas(H)=pas(i);
            sabias.missionid(H)=k;
            H=H+1;
        end
    end
end

%% plot bias
for k=1:3
figure,
latlim = [53 67];
lonlim = [10 31];
ax = usamap(latlim, lonlim);
hold on
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
setm(ax, 'FFaceColor', [1 1 1])
setm(gca,'FLineWidth',5,'Grid','on','FontSize',18,'fontweight','bold')
scatterm(sabias.lat(sabias.missionid==k),sabias.lon(sabias.missionid==k),30,sabias.bias(sabias.missionid==k),'filled','o');
c=colorbar; 
c.Label.String = 'Bias_S_A_-_H_D_M [cm]';
caxis([0 20]); %fix the bar: SA

set(c,'position',[.71 .2 .015 .7],'FontSize',18,'FontWeight','bold')


title(strcat(str(k,:),': ','Bias_S_A_-_H_D_M','Mean=',num2str(mean(sabias.bias(sabias.missionid==k)),2),' [cm]'),'FontSize',18,'FontWeight','bold');
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; ax.FontName='Times New Roman';
xLoc =1.3001e+05;
yLoc =7.0876e+06;
scaleruler('Units', 'km', 'RulerStyle', 'patches', ...
'XLoc', xLoc, 'YLoc', yLoc,"FontSize",10);

end

