function [base_state]= findStateHap(read_alignseq, quality_alignseq,mapq,ref_align,prim,tprob)

transitionmatrix=zeros(10,10);
transitionmatrix(1,1)=tprob(1,1);
transitionmatrix(1,2:4)=tprob(1,2)/3;
transitionmatrix(1,5)=tprob(1,3);
transitionmatrix(1,6:10)=tprob(1,4)/5;
transitionmatrix(2,1)=tprob(2,1);
transitionmatrix(2,2:4)=tprob(2,2)/3;
transitionmatrix(2,5)=tprob(2,3);
transitionmatrix(2,6:10)=tprob(2,4)/5;
transitionmatrix(3,:)=transitionmatrix(2,:);
transitionmatrix(4,:)=transitionmatrix(3,:);
transitionmatrix(5,1)=tprob(3,1);
transitionmatrix(5,2:4)=tprob(3,2)/3;
transitionmatrix(5,5)=tprob(3,3);
transitionmatrix(5,6:10)=tprob(3,4)/5;
transitionmatrix(6,1)=tprob(4,1);
transitionmatrix(6,2:4)=tprob(4,2)/3;
transitionmatrix(6,5)=tprob(4,3);
transitionmatrix(6,6:10)=tprob(4,4)/5;
transitionmatrix(7,:)=transitionmatrix(6,:);
transitionmatrix(8,:)=transitionmatrix(6,:);
transitionmatrix(9,:)=transitionmatrix(6,:);
transitionmatrix(10,:)=transitionmatrix(6,:);

base_state = zeros(1,size(ref_align,2));
%emission matrix

for i=1:size(ref_align,2)
    if size(read_alignseq(read_alignseq(:,i) ~= -1),1) < 5
        % -1 means not cover,  
        base_state (i) = -1;
        %hap_base(i) = -1;
    end
end
loc_sec1=find(diff(base_state)==1)+1;%beginning of each section
loc_sec2=find(diff(base_state)==-1);  %end of each section

for i=1:length(loc_sec1)
    ref_sec = ref_align(loc_sec1(i):loc_sec2(i));
    read_sec = read_alignseq(:,loc_sec1(i):loc_sec2(i));
    quality_sec = quality_alignseq(:,loc_sec1(i):loc_sec2(i));
    emissionmatrix = buildEmissionMatrixHap(ref_sec,read_sec,quality_sec,mapq,prim);
    emissionmatrix(isnan(emissionmatrix)) = -Inf;
    initialprobs=transitionmatrix(1,:).';
    mysequence=ref_sec;
    [Delta_sec, state_sec] = calculateDeltaHMM(mysequence, initialprobs, transitionmatrix, emissionmatrix);    
    base_state(loc_sec1(i):loc_sec2(i))=state_sec;
    %hap_base(loc_sec1(i):loc_sec2(i))=hap_sec;
end
end

    


