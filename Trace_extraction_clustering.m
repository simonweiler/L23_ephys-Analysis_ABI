%SW181215 Script
% Extract waveforms and features for electrophysiology similar to https://www.biorxiv.org/content/biorxiv/early/2018/07/17/368456.full.pdf
% In my data set, I only used 1 s long pulses 
% for calculating ap threshold etc.  I use https://senselab.med.yale.edu/SimToolDB/showTool.cshtml?tool=112112&file=%5Cpandora-1.3b%5Cpandora%5C%40spike_shape%5CcalcInitVmSekerliV2.m
%or check https://github.com/cengique/pandora-matlab

%Nested Functions:
%extract_ISI
%extract_Rheobase_traces
%ephys_parameter_bins
%MAKE SURE PANDORA TOOLBOX IS DOWNLOADED AND IN THE MAT PATH

%Allen Brain Institute Ephys Parameters:
%1.APVm: waveform of first AP; including 3ms after AP threshold 
%2.AP dV/dt: waveform; Time derivative of (1)
%3. ISI shape: waveform; Average of ISI voltage trajectories; aligned to the threshold of the initial AP and noralized in duration
%4: Subthreshold abs.: waveform; concanetad response to hyperpolarizing current
%steps--->>>> NOT POSSIBLE WITH MY DATA SET
%5: Subthreshold norm.: waveform; response to largest hyperpolarizing current
%step, aligned to baseline membrane potential and normalized by maximum
%voltage deflection 
%6: PSTH: AP counts in 50 ms bins; divided by bin width; from rheobase to
%rheobase +100 pA
%7: Inst.firing rate: Instantenous firing rate binned 20 ms;from rheobase to
%rheobase +100 pA
%8: Up/down: Upstroke/downstroke ratio; binned 20 ms;from rheobase to
%rheobase +100 pA
%9: AP peak: AP peak; binned 20 ms;from rheobase to
%rheobase +100 pA
%10: AP fast trough: AP fast trough;binned 20 ms;from rheobase to
%rheobase +100 pA---->>>>>> NO IDEA WHAT THEY MEAN
%11: AP thresh.: AP threshold; binned 20 ms;from rheobase to
%rheobase +100 pA 
%12: AP width: AP width;binned 20 ms;from rheobase to
%rheobase +100 pA 
%13: Inst. freq: Instantaneous firing rate normalized to maximum rate for
%each step;binned 20 ms;from rheobase to
%rheobase +100 pA 


%clean up
clear all;
close all;
%LOAD EPHYS STRUCTURE (containing the raw traces, and some event detection
pathName='R:\Share\Simon\Drago_Volker_Simon\Manuscript2018\Electrophysiology\Ephys structure';%NEEDS THE MAT FILE LOCATED THERE
list=dir([char(pathName) '\*.mat']);
load([char(pathName) '/' list.name],'-mat');
savefolder='R:\Share\Simon\Drago_Volker_Simon\Manuscript2018\Electrophysiology\Approach 3_ABI_feature_sPCA_GMM\out'
%--------------------------
%IMPORTANT FLAGS
disp1_2=0;%plot or not parameters 1 & 2
par1_2=0;%calculate/extract parameters 1 & 2
disp3=1;%same logic 
par3=1;
disp5=0;
par5=0;
par6_13=1;
disp6_13=1;
%%%%%%%%%%%%%%%%%

L23_ephys={};%empty structure for saving variables

%Preprocess data: remove nonsense
%in total 161 cells: 
%1-13 cells: 20 pA and 40 pA step missing; after 150 pA in steps of 50 pA:
%current=[10 30 50 60 70 80 90 100 110 120 130 140 150 200 250 300 350]
%(17 traces in total)
%%%%%%%%%%%%%%%
%14-90 cells:20 pA 40 pA step included; after 150 pA in steps of 50 pA 
%current=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 200 250 300 350]
%(19 traces in total)
%%%%%%%%%%%%%%%
%91-161 cells
%%current=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 180 200 220 240 260 280 300 400]
%(24 traces)


all=ephys.all;%all is extracted features
all=ephys.all';
m=[4;9;10;13;16;17];%remove unusable data (containing NaN or other issues);
all(:,m)=[];
vhigh=find(all(:,13)>-55);%remove cells more depolarized than -55 mV (quality criterion)
all(vhigh,:)=[];%1-13;14-86;87-144
%---------------------------
%spike time for each trial for each cell 
event=ephys.spike_event;
event(:,vhigh)=[];
%Rheobase in pA for each cell 
rheobase_v=ephys.Rheobase;
rheobase_v(:,vhigh)=[];
%raw ephys traces
traces=ephys.spike_traces; 
traces(:,:,vhigh)=[];
traces_all=ephys.ephystraces;
traces_all(:,:,vhigh)=[];
% 144 cells in total 
%% 
%%Extracting Parameters
%%PARAMETERS: 1, 2 %1.APVm; %2.AP dV/dt
if par1_2==1;
%AP waveform parameter 1 and 2 
for nr=1:length(event)
ap(:,nr)=find(~cellfun(@isempty,event(:,nr)),1);%find indx where Rheobase occurs
tr=cell2mat(traces(:,:,nr));
startp(:,nr)=find(tr(:,ap(:,nr),1)>=all(nr,3),1);%find ephys trace where first AP occurs
end
nr=[];
for i=1:length(event)
APwave=cell2mat(traces(:,:,i));
APwaveform(:,i)=APwave(startp(:,i):startp(:,i)+300,ap(:,i));%%%%%%%  PARAMETER 1 %%%%%%%%
APwaveform_der(:,i)=diff(APwave(startp(:,i):startp(:,i)+300,ap(:,i)));%%%%%%%   PARAMETER 2    %%%%%%%%% 
end
i=[];
ap=[];
tr=[];
%PLOT
if disp1_2==1;
f1=figure('Name','APVm');
set(gcf, 'Position', [200, 0, 1500, 1000]);
f2=figure('Name','AP dV/dt');
set(gcf, 'Position', [200, 0, 1500, 1000]);
figure(f1);
for i=1:length(event)
subplot(12,12,i);
imagesc(APwaveform(:,i)');
hold on;
%xlabel('Time');
set(gca,'Yticklabel',[]);
%colorbar;
caxis([-72 78]);
end
i=[];
figure(f2);
for i=1:length(event)
subplot(12,12,i);
imagesc(APwaveform_der(:,i)');
hold on;
%xlabel('Time');
set(gca,'Yticklabel',[]);
%colorbar;
caxis([-8 24]);
end
i=[];
end 
%Put in structure 
L23_ephys.data{1,1}='APVm';
L23_ephys.data{1,2}=APwaveform;%PARAMETER 1
L23_ephys.data{2,1}='AP dV/dt';
L23_ephys.data{2,2}=APwaveform_der;%PARAMETER 2
end
%% %%PARAMETER: 3 ISI_shape
if par3==1;
%call function extract_ISI
[ISI_final]=extract_ISI(event, traces);
%PLOT
if disp3==1
f4=figure('Name','ISI shape');set(gcf, 'Position', [200, 0, 1500, 1000]);
title('ISI shape')
subplot(1,2,1);
plot(ISI_final);
xlabel('Time');
ylabel('delta Voltage (mV)');
axis square;
title('ISI shape');
hold on;
subplot(1,2,2);
imagesc(ISI_final');
xlabel('Time');
ylabel('Cell');
colorbar;
axis square;
end
L23_ephys.data{3,1}='ISI shape';
L23_ephys.data{3,2}=ISI_final;%PARAMETER 3
end
%% 
%%PARAMETERS: 5; 
if par5==1;
for i=1:size(traces_all,3)
[a,b,c]=find(traces_all(:,:,i)==min(min(traces_all(:,:,i))));
sag=traces_all(:,b,i);
sag_s=sag-mean(sag(1:1000));%Baseline subtracted  (starts at 0)
sag_n=sag_s/-min(sag_s);%normalized to peak
if ~isempty(find(sag_n(1500:10000)>-0.8));
    sag=[];
    sag_s=[];
    sag_n=[];
    sag=traces_all(:,b+1,i);
    sag_s=sag-mean(sag(1:1000));%Baseline subtracted  (starts at 0)
    sag_n=sag_s/-min(sag_s);%normalized to peak
end
sag_traces(:,i)=sag_n(1:12000,:);;%%%%%%%   PARAMETER 5    %%%%%%%%% 
end
i=[];
%PLOT
if disp5==1
f3=figure('Name','Subthreshold norm');
for i=1:size(traces_all,3)
hold on;
plot(sag_traces(:,i));
%imagesc(sag_traces');colorbar;
xlabel('Time');
ylabel('Norm. amplitude');
end
end
i=[];
%Put in structure 
L23_ephys.data{4,1}='Subthreshold norm';
L23_ephys.data{4,2}=sag_traces;%PARAMETER 5
end
%% PARAMETERS: 6, 7, 11, 12, 13
if par6_13==1;
%call function extract_Rheobase_traces    
[traces_Rheo]=extract_Rheobase_traces(rheobase_v, traces, event);
%PLOT
if disp6_13==1
f6=figure;set(gcf, 'Position', [200, 0, 1500, 1000]);
for i=1:size(traces_Rheo,3)
subplot(12,12,i)
hold on;
plot(traces_Rheo(:,:,i));
end
end
% EXTRACT PARAMETERS, time in ms and voltage in mV
%call function ephys_parameter_bins
[AP_thresh, AP_peak, AP_width, AP_ifreq, AP_ifreq_n, AP_PSTH, AP_upstroke, AP_downstroke]= ephys_parameter_bins(traces_Rheo,size(traces_Rheo,3));
AP_up_down_ratio=AP_upstroke./AP_downstroke;
matr={AP_thresh, AP_peak, AP_width, AP_ifreq, AP_ifreq_n, AP_up_down_ratio};
%replicate values wehn lees that 6 sweeps
for k=1:6
for i=1:length(AP_thresh);
    temp=matr{k};
    inter=temp(:,:,i);
    idx_col=find(isnan(temp(1,:,i)));
    if isempty(idx_col)
      temp_r(:,:,i)=inter;
    else
    idx_rep=idx_col(1)-1;
    col_rep=temp(:,idx_rep,i);
    rep=repmat(col_rep,1,length(idx_col));
    inter(:,idx_col)=rep;
    temp_r(:,:,i)=inter;
    
    end
    
end
temp_reshap=reshape(temp_r,[],length(AP_thresh),1);
comb(:,:,k)=temp_reshap;
end
  
AP_PSTH_r=reshape(AP_PSTH,[],144,1);

if disp6_13==1
f7=figure;set(gcf, 'Position', [200, 0, 1500, 1000]);
for i=1:size(traces_Rheo,3)
subplot(12,12,i)
hold on;
plot(AP_thresh_r(:,:,i));
end
end
L23_ephys.data{5,1}='PSTH';
L23_ephys.data{5,2}=AP_PSTH_r;%PARAMETER 6
L23_ephys.data{6,1}='Inst.firing rate';
L23_ephys.data{6,2}=comb(:,:,4);%PARAMETER 7
L23_ephys.data{7,1}='Up/down';
L23_ephys.data{7,2}=comb(:,:,6);%PARAMETER 8
L23_ephys.data{8,1}='AP_peak';
L23_ephys.data{8,2}=comb(:,:,2);%PARAMETER 9
L23_ephys.data{9,1}='AP_thresh';
L23_ephys.data{9,2}=comb(:,:,1);%PARAMETER 11
L23_ephys.data{10,1}='AP_width';
L23_ephys.data{10,2}=comb(:,:,3);%PARAMETER 12
L23_ephys.data{11,1}='Norm.Inst.firing rate';
L23_ephys.data{11,2}=comb(:,:,5);%PARAMETER 13
%AP trough is missing -> DISCUSS

% SAVE in savefolder directory   
cd(savefolder);
FileName=['Data_SW_L23ephys',datestr(now, 'hh-dd-mmm-yyyy')];
save(FileName,'-struct','L23_ephys');
end

