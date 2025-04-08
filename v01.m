% make sa data cell based on said
for q=1:3
    SA{q}=sa(sa.said==q,:)
end
%% 

for k=1:3
    meanSA{k}=table();
    H=1;

    id=unique(SA{k}.id);
    for i=1:length(id)
        temp=SA{k}(SA{k}.id==id(i),:);
        for j=1:length(Lat1)
            for l=1:length(Lon1)
                samean=mean(temp.DTsa(temp.Lat>=Lat1(j)-deltalat/2&temp.Lat<Lat1(j)+deltalat/2&...
                    temp.Lon>=Lon1(l)-deltalon/2&temp.Lon<Lon1(l)+deltalon/2),'omitnan');

                if ~isnan(samean)
                    hdmmean=mean(temp.DThdm(temp.Lat>=Lat1(j)-deltalat/2&temp.Lat<Lat1(j)+deltalat/2&...
                        temp.Lon>=Lon1(l)-deltalon/2&temp.Lon<Lon1(l)+deltalon/2),'omitnan');

                    meanSA{k}.sa(H,1)=samean;
                    meanSA{k}.hdm(H,1)=hdmmean;
                    meanSA{k}.lat(H,1)=Lat1(j);
                    meanSA{k}.lon(H,1)=Lon1(l);
                    meanSA{k}.pas(H,1)=temp.Pass(1);
                    meanSA{k}.cycle(H,1)=temp.Cycle(1);
                    meanSA{k}.time(H,1)=temp.Time(1);
                    H=H+1;
                end
                clearvars hdmmean samean
            end
        end
        clear temp
    end
end
clearvars i j k l id temp