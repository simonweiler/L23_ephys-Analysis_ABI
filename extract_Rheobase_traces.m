function [traces_Rheo]=extract_Rheobase_traces(Rheobase, traces, event)%SW181215
%Output: 
%traces_Rheo: sorted Rheobase +100 pA in 20 pA steps traces per cell
%Input: 
%Rheobase: Rheobase in pA for each cell 
%traces: spike traces
%event: spike detected event with time 

%get trace at Rheobase to Rheobase+100 pA in 20 pA steps (so maximum 6
%traces per cell)-> PROBLEM: cells have different cuttent step size etc. Some cells have not all traces (max. 6) if not 6 then NaN as replacement 
%split cells: 
%three cell groups: c1, c2, c3
%current cell group 1=[10 30 50 60 70 80 90 100 110 120 130 140 150 200 250 300 350]
%current cell group 2=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 200 250 300 350]
%current cell group 3=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 180 200 220 240 260 280 300 400];
%three groups of cells based on their current step
c1=1:13;
c2=14:86;
c3=87:144;
%create matrix with zeros 
traces_Rheo=zeros([15000 6 144]);
%find the index for the Rheobase starting pulse 
for nr=1:length(event)
ap(:,nr)=find(~cellfun(@isempty,event(:,nr)),1);
end
%used current steps (in pA)
cu_step=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 180 200 220 240 260 280 300 400];
cu_step2=[10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 200 250 300 350 400 0 0 0 0];
%c1 idx 
for i=1:length(cu_step2)
    idx_c1{:,i}=find(Rheobase(c1)==cu_step2(i))
end
%c2 idx 
for i=1:length(cu_step)
    idx_c2{:,i}=find(Rheobase(c2)==cu_step2(i))
end
%c3 idx
for i=1:length(cu_step)
    idx_c3{:,i}=find(Rheobase(c3)==cu_step(i))
end
i=[];
%%%% 6 TRACES %%%%%
%starting from Rheobase 50 for c1 and c2; from Rheobase 70 for c3; then move in steps of 20 pA to 150 pA
idx6=[cell2mat(idx_c1(1:5)) cell2mat(idx_c2(1:5))+13 cell2mat(idx_c3(1:7))+86];
%%%% 5 TRACES %%%%%
idx5=[cell2mat(idx_c1(6:7)) cell2mat(idx_c2(6:7))+13];
%use idx and find traces
%%%% 4 TRACES %%%%%
idx4=[cell2mat(idx_c1(8:9)) cell2mat(idx_c2(8:9))+13];
%%%% 3 TRACES %%%%%
idx3=[cell2mat(idx_c1(10:11)) cell2mat(idx_c2(10:11))+13];
%%%% 2 TRACES %%%%%
idx2=[cell2mat(idx_c1(12:13)) cell2mat(idx_c2(12:13))+13];
%%%% 1 TRACE %%%%%
idx1=[cell2mat(idx_c1(14:24)) cell2mat(idx_c2(14:24))+13];
idx1_6={idx6 idx5 idx4 idx3 idx2 idx1};
%use idx and find traces for c1 and c2-> 6-1 traces
for i=1:6;
idx=idx1_6{i}    
a=ap(:,idx);
a_length(i)=length(a);
end
idx=idx1_6{1}
a=ap(:,idx);
for nr=1:a_length(1)
    tr=cell2mat(traces(:,:,idx(nr)));
    traces_Rheo_a6(:,:,nr)=tr(:,a(:,nr):2:a(:,nr)+11);
end
idx=idx1_6{2}
a=ap(:,idx);
for nr=1:a_length(2)
    tr=cell2mat(traces(:,:,idx(nr)));
    traces_Rheo_a5(:,:,nr)=tr(:,a(:,nr):2:a(:,nr)+9);  
end
idx=idx1_6{3}
a=ap(:,idx);
for nr=1:a_length(3)
    tr=cell2mat(traces(:,:,idx(nr)));
    traces_Rheo_a4(:,:,nr)=tr(:,a(:,nr):2:a(:,nr)+7);
end
idx=idx1_6{4}
a=ap(:,idx);
for nr=1:a_length(4)
    tr=cell2mat(traces(:,:,idx(nr)));
    traces_Rheo_a3(:,:,nr)=tr(:,a(:,nr):2:a(:,nr)+5);
end
idx=idx1_6{5}
a=ap(:,idx);
for nr=1:a_length(5)
    tr=cell2mat(traces(:,:,idx(nr)));
    traces_Rheo_a2(:,:,nr)=tr(:,a(:,nr):2:a(:,nr)+3);
end 
idx=idx1_6{6}
a=ap(:,idx);
for nr=1:a_length(6)
    tr=cell2mat(traces(:,:,idx(nr)));
    traces_Rheo_a1(:,:,nr)=tr(:,a(:,nr));
end 
%-->101 cells, 86 c1 and c2, 15 c3

%for c3 find further idx with Rheobase 80 or more   
for i=8:24
idx_c3_6{:,i}=cell2mat(idx_c3(i))+86%-> 43 cells
end
%6 traces
%find traces for Rheobase 80 
idx6_8=idx_c3_6{:,8};
a6_8=ap(:,idx6_8);
for nr=1:length(a6_8)
    tr=cell2mat(traces(:,:,idx6_8(nr)));
    traces_Rheo_a6_8(:,:,nr)=tr(:,[a6_8(:,nr):2:a6_8(:,nr)+9 a6_8(:,nr)+10]);
end
nr=[];
tr=[];
%find traces for Rheobase 90 
idx6_9=idx_c3_6{:,9};
a6_9=ap(:,idx6_9);
for nr=1:length(a6_9)
    tr=cell2mat(traces(:,:,idx6_9(nr)));
    traces_Rheo_a6_9(:,:,nr)=tr(:,[a6_9(:,nr):2:a6_9(:,nr)+8 a6_9(:,nr)+9]);
end
nr=[];
tr=[];
%find traces for Rheobase 100
idx6_10=idx_c3_6{:,10};
a6_10=ap(:,idx6_10);
for nr=1:length(a6_10)
    tr=cell2mat(traces(:,:,idx6_10(nr)));
    traces_Rheo_a6_10(:,:,nr)=tr(:,[a6_10(:,nr):2:a6_10(:,nr)+7 a6_10(:,nr)+8:a6_10(:,nr)+9]);
end
nr=[];
tr=[];

%find traces for Rheobase 110
idx6_11=idx_c3_6{:,11};
a6_11=ap(:,idx6_11);
for nr=1:length(a6_11)
    tr=cell2mat(traces(:,:,idx6_11(nr)));
    traces_Rheo_a6_11(:,:,nr)=tr(:,[a6_11(:,nr):2:a6_11(:,nr)+5 a6_11(:,nr)+7:a6_11(:,nr)+9]);
end
nr=[];
tr=[];
%find traces for Rheobase 120
idx6_12=idx_c3_6{:,12};
a6_12=ap(:,idx6_12);
for nr=1:length(a6_12)
    tr=cell2mat(traces(:,:,idx6_12(nr)));
    traces_Rheo_a6_12(:,:,nr)=tr(:,[a6_12(:,nr):2:a6_12(:,nr)+5 a6_12(:,nr)+6:a6_12(:,nr)+8]);
end
nr=[];
tr=[];
%find traces for Rheobase 130
idx6_13=idx_c3_6{:,13};
a6_13=ap(:,idx6_13);
for nr=1:length(a6_13)
    tr=cell2mat(traces(:,:,idx6_13(nr)));
    traces_Rheo_a6_13(:,:,nr)=tr(:,[a6_13(:,nr):2:a6_13(:,nr)+5 a6_13(:,nr)+6:a6_13(:,nr)+8]);
end
nr=[];
tr=[];
%find traces for Rheobase 140
idx6_14=idx_c3_6{:,14};
a6_14=ap(:,idx6_14);
for nr=1:length(a6_14)
    tr=cell2mat(traces(:,:,idx6_14(nr)));
    traces_Rheo_a6_14(:,:,nr)=tr(:,[a6_14(:,nr):2:a6_14(:,nr)+3 a6_14(:,nr)+4:a6_14(:,nr)+7]);
end
nr=[];
tr=[];
idx6_15=[idx_c3_6{15:17}];
a6_15=ap(:,idx6_15);
for nr=1:length(a6_15)
    tr=cell2mat(traces(:,:,idx6_15(nr)));
    traces_Rheo_a6_15(:,:,nr)=tr(:,[a6_15(:,nr):1:a6_15(:,nr)+5]);
end
nr=[];
tr=[];
%5 traces c3
idx5_2=[idx_c3_6{18}];
a5_2=ap(:,idx5_2);
for nr=1:length(a5_2)
    tr=cell2mat(traces(:,:,idx5_2(nr)));
    traces_Rheo_a5_2(:,:,nr)=tr(:,[a5_2(:,nr):1:a5_2(:,nr)+4]);
end
nr=[];
tr=[];
%4 traces c3->empty
%3 traces c3->empty
%2 traces c3-> empty
%1 trace c3->empty

%combine
traces_Rheo_6=cat(3,traces_Rheo_a6,traces_Rheo_a6_8,traces_Rheo_a6_9,traces_Rheo_a6_10, traces_Rheo_a6_11, traces_Rheo_a6_12, traces_Rheo_a6_13, traces_Rheo_a6_14, traces_Rheo_a6_15);
traces_Rheo_5=cat(3,traces_Rheo_a5,traces_Rheo_a5_2);
id6=[idx6 idx6_8 idx6_9 idx6_10 idx6_11 idx6_12 idx6_13 idx6_14 idx6_15];
id5=[idx5 idx5_2];
traces_Rheo_5=[traces_Rheo_5 NaN([15000 1 length(id5)])];
traces_Rheo_4=[traces_Rheo_a4 NaN([15000 2 length(idx4)])];
traces_Rheo_3=[traces_Rheo_a3 NaN([15000 3 length(idx3)])];
traces_Rheo_2=[traces_Rheo_a2 NaN([15000 4 length(idx2)])];
traces_Rheo_1=[traces_Rheo_a1 NaN([15000 5 length(idx1)])];
for nr=1:length(id6)
    traces_Rheo(:,:,id6(nr))=traces_Rheo_6(:,:,nr);
end
nr=[];
for nr=1:length(id5)
    traces_Rheo(:,:,id5(nr))=traces_Rheo_5(:,:,nr);
end
nr=[];
for nr=1:length(idx4)
    traces_Rheo(:,:,idx4(nr))=traces_Rheo_3(:,:,nr);
end
nr=[];
for nr=1:length(idx3)
    traces_Rheo(:,:,idx3(nr))=traces_Rheo_3(:,:,nr);
end
nr=[];
for nr=1:length(idx2)
    traces_Rheo(:,:,idx2(nr))=traces_Rheo_2(:,:,nr);
end
nr=[];
for nr=1:length(idx1)
    traces_Rheo(:,:,idx1(nr))=traces_Rheo_1(:,:,nr);
end

end