function [ISI]=extract_ISI(event, traces)
%Output: ISI->waveform; Average of ISI voltage trajectories for each cell; aligned to the threshold of the initial AP and noralized in duration
%Input
%spikes: number of spikes per trace per cell 
%traces: time based spike traces
spikes=cellfun(@length,event);
spike_count=max(cellfun(@length,event));
AP5=find(spike_count>=5);
AP5_n=find(spike_count<5);
dt = 1e-4;
dy = 1e-3;
lAP={AP5_n' AP5'};

for jj=1:2
    AP=cell2mat(lAP(jj))';
for i=1:length(AP)
tr5=cell2mat(traces(:,:,AP(i)));
if jj==1;
nr=find((spikes(:,AP(i)))>=max(spikes(:,AP(i))));
nr_all(:,i)=length(find((spikes(:,AP(i)))>=max(spikes(:,AP(i)))));
%nr=nr(end);
else
nr=[];
nr=find(spikes(:,AP(i)));
nr_all(:,i)=1;
nr=nr(1);
end
tr5_all{:,i}=tr5(:,nr);
for t=1:size(cell2mat(tr5_all(i)),2)
trace_active=cell2mat(tr5_all(i));
props = struct('spike_finder', 2, 'threshold', 0);
traces_analysis(:,t)=trace(trace_active(:,t),dt, dy, 'Analysis', props);
alltrace_info(:,t) = getProfileAllSpikes(traces_analysis(:,t));
try
parameters{:,t,i}=alltrace_info(1,t).spikes_db.data;
catch
parameters{:,t,i}=NaN;
end

for tt=1:size(parameters{:,t,i},1)-1   
if length(alltrace_info(1,t).spikes(1,:).times)>=2;
temp=parameters{:,t,i};
time=alltrace_info(1,t).spikes.times;
try
loc=find(trace_active(time(tt):time(tt+1),t)<=temp(tt,1));
ftime=loc(1)+time(tt);
end
try
loc2=find(trace_active(time(tt):time(tt+1),t)<=temp(tt+1,4));
ftime2=loc2(end)+time(tt);
end 
else
temp=parameters{:,t,i};
time=alltrace_info(1,t).spikes.times;
loc=find(trace_active(time(tt):time(tt)+1000,t)<=temp(tt,1));
ftime=loc(1)+time(tt);
loc2=find(trace_active(time(tt):time(tt)+1000,t)<=temp(tt,4));
ftime2=loc2(end)+1000;

end
ISI{:,tt,t,i}=trace_active(ftime:ftime2,t)-temp(tt,4);
tempv=trace_active(ftime:ftime2,t)-temp(tt,4);
tempISI(:,tt)=length(trace_active(ftime:ftime2,t)-temp(tt,4)); 
end
end
end

for i=1:length(AP)
 for t=1:size(cell2mat(tr5_all(i)),2)
   for tt=1:size(parameters{:,t,i},1)
      try
        oldx=[1:length(cell2mat(ISI(:,tt,t,i)))];
        oldy=cell2mat(ISI(:,tt,t,i));
           b=max(cellfun('length',ISI));
           newx=[1:b(:,:,t,i)];
           xval=linspace(1,b(:,:,t,i),length(cell2mat(ISI(:,tt,t,i))));
           newY = interp1(xval, oldy, newx,'linear');   
       end
        try
            oldx_d=newx;
            oldy_d=newY;
            newx_d=[1:b(:,:,t,i)/100:b(:,:,t,i)];
             newY_d = interp1(oldx_d, oldy_d, newx_d);
        end      
         ISI_standard{:,tt,t,i}=newY_d;      
            
end 
end
end
i=[];
for i=1:length(AP)
    
    temp_t=ISI_standard(:,:,nr_all(i),i);
    ISI_ave2(:,i)=nanmean(cat(3,temp_t{:}),3);
    if isnan(ISI_ave2(:,i))==1
     temp_t=ISI_standard(:,:,nr_all(i)-1,i);
    ISI_ave2(:,i)=nanmean(cat(3,temp_t{:}),3);   
    else
    ISI_ave2(:,i)=nanmean(cat(3,temp_t{:}),3);
    end
    ISI_ave2_smooth(:,i)=smooth(ISI_ave2(:,i),0.1,'rloess');
    
end
ISI_t{:,jj}=ISI_ave2_smooth;
end
ISI=zeros(100,length(AP)+length(AP5_n));
ISI(:,AP5_n)=cell2mat(ISI_t(1));
ISI(:,AP5)=cell2mat(ISI_t(2));
end