function [thresh, peak, width, ifreq, ifreq_n, PSTH, upstroke, downstroke]= ephys_parameter_bins(traces_Rheo,nr_cells) 
%use the Rheobase to +100 pA traces and extract spike times and theshold of
%spikes in a bin of 20 ms (51 bins); if there is no spike in a bin interpolate/extrapolate from nearest
%if sweep is missing NaN

%OUTPUT: all binned values (either 50 or 20 ms bin per trace=
%thresh: AP threshold (in mV)
%peak: AP peak amplitude (in mV)
%width: AP width (in ms)
%ifreq: inst. firing frequency (in Hz)
%ifreq_n= Instantaneous firing rate normalized to maximum rate per step
%PSTH: PSTH
%upstroke: AP max dv/dt
%downstroke: AP min dv/dt

%INPUT
%traces_Rheo: sorted Rheobase +100 pA in 20 pA steps traces per cell
%nr_cells: number of cells to analyze 


%for calculating ap threshold etc.  I use https://senselab.med.yale.edu/SimToolDB/showTool.cshtml?tool=112112&file=%5Cpandora-1.3b%5Cpandora%5C%40spike_shape%5CcalcInitVmSekerliV2.m
%MAKE SURE PANDORA TOOLBOX IS DOWNLOADED AND IN THE MAT PATH
%entries 
dt = 1e-4;
dy = 1e-3;
for p=1:nr_cells%cell number
trace_active=traces_Rheo(:,:,p);
props = struct('spike_finder', 2, 'threshold', 0);
for i=1:size(trace_active,2)%sweep number
trace_curr=trace_active(:,i);
traces_analysis=trace(trace_active(:,i),dt, dy, 'Analysis', props);%function from PANDORA TOOLBOX
alltrace_info = getProfileAllSpikes(traces_analysis);%function from PANDORA TOOLBOX
for k=1:size(alltrace_info(1,1).spikes_db(:,:).data,1);%function from PANDORA TOOLBOX; how many spikes per sweep
parameters=alltrace_info(1,1).spikes_db(k,:).data;%function from PANDORA TOOLBOX
%Parameters are:
% %%%%[ 1]    'MinVm'        
%     [ 2]    'PeakVm'       
%     [ 3]    'InitVm'       
%     [ 4]    'InitVmBySlope'
%     [ 5]    'MaxVmSlope'   
%     [ 6]    'HalfVm'       
%     [ 7]    'Amplitude'    
%     [ 8]    'MaxAHP'       
%     [ 9]    'DAHPMag'      
%     [10]    'InitTime'     
%     [11]    'RiseTime'     
%     [12]    'FallTime'     
%     [13]    'MinTime'      
%     [14]    'BaseWidth'    
%     [15]    'HalfWidth'    
%     [16]    'FixVWidth'    
%     [17]    'Index'        
%     [18]    'Time'     
time_ap(:,k)=parameters(18);%extract time of each spike
thresold_ap(:,k)=parameters(3);%spike threshold, parameter 3
peak_ap(:,k)=parameters(2);%spike peak, parameter 3
hwidth_ap(:,k)=parameters(15);%half width, parameter 15
ISI_ap(:,k)=1./diff([0 time_ap(:,k)]);%1/ISI
if any(isnan(thresold_ap)) % if this is true, run the 'continue'
    continue
end
slope_up(:,k)=max(diff(trace_curr(round(time_ap(:,k)*10-20):round(time_ap(:,k)*10+30))))*10;
slope_down(:,k)=min(diff(trace_curr(round(time_ap(:,k)*10-20):round(time_ap(:,k)*10+30))))*10;
end

if size(alltrace_info(1,1).spikes_db(:,:).data,1)>1
try
vi(:,i)=interp1(time_ap,thresold_ap,[100:20:1100],'nearest','extrap');%interpolate in 51 bins
vp(:,i)=interp1(time_ap,peak_ap,[100:20:1100],'nearest','extrap');
vw(:,i)=interp1(time_ap,hwidth_ap,[100:20:1100],'nearest','extrap');
vf(:,i)=interp1(time_ap,ISI_ap,[100:20:1100],'nearest','extrap');
vsu(:,i)=interp1(time_ap,slope_up,[100:20:1100],'nearest','extrap');
vsd(:,i)=interp1(time_ap,slope_down,[100:20:1100],'nearest','extrap');
catch
vi(:,i)=NaN(1,51);%interpolate in 51 bins
vp(:,i)=NaN(1,51);
vw(:,i)=NaN(1,51);
vf(:,i)=NaN(1,51);
vsu(:,i)=NaN(1,51);
vsd(:,i)=NaN(1,51);
end
else
try
vi(:,i)=repmat(thresold_ap(1),1,51);%simply replicate if there is only one spike in entire sweep
vp(:,i)=repmat(peak_ap(1),1,51);
vw(:,i)=repmat(hwidth_ap(1),1,51);
vf(:,i)=repmat(ISI_ap(1),1,51);
vsu(:,i)=repmat(slope_up(1),1,51);
vsd(:,i)=repmat(slope_down(1),1,51);
catch
vi(:,i)=NaN(1,51);%interpolate in 51 bins
vp(:,i)=NaN(1,51);
vw(:,i)=NaN(1,51);
vf(:,i)=NaN(1,51);
vsu(:,i)=NaN(1,51);
vsd(:,i)=NaN(1,51);
end
end
vf2(:,i)=histcounts(time_ap-100,20)./50;
vf_n(:,i)=vf(:,i)./max(vf(:,i));
end
if size(vi,1)<6 %simply extrapolate values when less than 6 sweeps; not happy with that
%      for t=(size(vi,1)+1):6
         for zz=1:51
         vi_extr(:,zz)=interp1([1:size(vi,1)],vi(:,zz),[size(vi,1)+1:6],'linear','extrap');
         vp_extr(:,zz)=interp1([1:size(vp,1)],vp(:,zz),[size(vp,1)+1:6],'linear','extrap');
         vw_extr(:,zz)=interp1([1:size(vw,1)],vw(:,zz),[size(vw,1)+1:6],'linear','extrap');
         vf_extr(:,zz)=interp1([1:size(vf,1)],vf(:,zz),[size(vf,1)+1:6],'linear','extrap');
         vf_n_extr(:,zz)=interp1([1:size(vf_n,1)],vf_n(:,zz),[size(vf_n,1)+1:6],'linear','extrap');
         vsu_extr(:,zz)=interp1([1:size(vsu,1)],vsu(:,zz),[size(vsu,1)+1:6],'nearest','extrap');
         vsd_extr(:,zz)=interp1([1:size(vsd,1)],vsd(:,zz),[size(vsd,1)+1:6],'nearest','extrap');
         end
         for zz=1:20
         vf2_extr(:,zz)=interp1([1:size(vf2,1)],vf2(:,zz),[size(vf2,1)+1:6],'linear','extrap');
         end
         vi=[vi;vi_extr];
         vp=[vp;vp_extr];
         vw=[vw;vw_extr];
         vf=[vf;vf_extr];
         vf2=[ vf2; vf2_extr];
         vf_n=[vf_n;vf_n_extr];
         vsu=[vsu;vsu_extr];
         vsd=[vsd;vsd_extr];
         
%     vi(t,:)=vi(end,:);
%     vp(t,:)=vp(end,:); 
%     vw(t,:)=vw(end,:);
%      vf(t,:)=vf(end,:);
%      vf2(t,:)=vf2(end,:);
%      vf_n(t,:)=vf_n(end,:);
%      end
else
     vi=vi;
     vp=vp;
     vw=vw;
     vf=vf;
     vf2=vf2;
     vf_n=vf_n;
     vsu=vsu;
     vsd=vsd;
end
 
time_ap=[];
thresold_ap=[];
peak_ap=[];
hwidth_ap=[];
ISI_ap=[];
slope_up=[];
slope_down=[];
thresh{:,p}=vi;%%%%%%%   PARAMETER 11   %%%%%%%%%
peak{:,p}=vp;%%%%%%%   PARAMETER 9    %%%%%%%%%
width{:,p}=vw;%%%%%%%   PARAMETER 12    %%%%%%%%%
ifreq{:,p}=vf;%%%%%%%   PARAMETER 7    %%%%%%%%%
PSTH{:,p}=vf2;%%%%%%%   PARAMETER 6    %%%%%%%%%
ifreq_n{:,p}=vf_n;%%%%%%%   PARAMETER 13    %%%%%%%%%
upstroke{:,p}=vsu;
downstroke{:,p}=vsd;

vi_extr=[];
vp_extr=[];
vw_extr=[];
vf_extr=[];
vf2_extr=[];
vf_n_extr=[];
vsu_extr=[];
vsd_extr=[];
vi=[];
vp=[];
vw=[];
vf=[];
vf2=[];
vf_n=[];
vsu=[];
vsd=[];
clearvars parameters;
end
thresh=reshape(cell2mat(thresh),[51,6,p]);
peak=reshape(cell2mat(peak),[51,6,p]);
width=reshape(cell2mat(width),[51,6,p]);
ifreq=reshape(cell2mat(ifreq),[51,6,p]);
ifreq_n=reshape(cell2mat(ifreq_n),[51,6,p]);
PSTH=reshape(cell2mat(PSTH),[20,6,p]);
upstroke=reshape(cell2mat(upstroke),[51,6,p]);
downstroke=reshape(cell2mat(downstroke),[51,6,p]);


end