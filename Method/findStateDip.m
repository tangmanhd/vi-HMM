function [base_state]= findStateDip(read_alignseq, quality_alignseq,mapq,ref_align,prim,tprob,hetrate)
%tranistion matrix
m=30;
transitionmatrix=zeros(m,m);
transitionmatrix(1,:)= buildTrans(tprob(1,:),hetrate);
transitionmatrix(2,:)= buildTrans(tprob(2,:),hetrate);
transitionmatrix(3,:) = transitionmatrix(2,:);
transitionmatrix(4,:) = transitionmatrix(2,:);
transitionmatrix(5,:) = transitionmatrix(2,:);
transitionmatrix(6,:) = transitionmatrix(2,:);
transitionmatrix(7,:) = transitionmatrix(2,:);
transitionmatrix(9,:) = transitionmatrix(2,:);
transitionmatrix(10,:) = transitionmatrix(2,:);
transitionmatrix(12,:) = transitionmatrix(2,:);
transitionmatrix(15,:) = buildTrans(tprob(3,:),hetrate);
transitionmatrix(8,:) = (transitionmatrix(2,:)+transitionmatrix(15,:))/2;
transitionmatrix(11,:) = transitionmatrix(8,:);
transitionmatrix(13,:) = transitionmatrix(8,:);
transitionmatrix(14,:) = transitionmatrix(8,:);
transitionmatrix(16,:) = buildTrans(tprob(4,:),hetrate);
transitionmatrix(17,:) = transitionmatrix(16,:);
transitionmatrix(18,:) = transitionmatrix(16,:);
transitionmatrix(19,:) = transitionmatrix(16,:);
transitionmatrix(20,:) = transitionmatrix(16,:);
transitionmatrix(21,:) = transitionmatrix(16,:);
transitionmatrix(22,:) = transitionmatrix(16,:);
transitionmatrix(23,:) = transitionmatrix(16,:);
transitionmatrix(24,:) = transitionmatrix(16,:);
transitionmatrix(25,:) = transitionmatrix(16,:);
transitionmatrix(26,:) = transitionmatrix(16,:);
transitionmatrix(27,:) = transitionmatrix(16,:);
transitionmatrix(28,:) = transitionmatrix(16,:);
transitionmatrix(29,:) = transitionmatrix(16,:);
transitionmatrix(30,30) = transitionmatrix(29,30);
transitionmatrix(30,1:29) = (1-transitionmatrix(30,30))/29;

base_state = zeros(1,size(ref_align,2));
%hap_base = zeros(1,size(alignmat,2));% store haplotype
%emission matrix

for i=1:size(ref_align,2)
    if size(read_alignseq(read_alignseq(:,i) ~= -1),1) < 5 || ref_align(i) == 7
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
    emissionmatrix = buildEmissionMatrixDip(ref_sec,read_sec,quality_sec,mapq,prim);
    emissionmatrix(isnan(emissionmatrix)) = -Inf;
    initialprobs=transitionmatrix(1,:).';
    mysequence=ref_sec;
    [Delta_sec, state_sec] = calculateDeltaHMM(mysequence, initialprobs, transitionmatrix, emissionmatrix);
    base_state(loc_sec1(i):loc_sec2(i))=state_sec;
end
end

    
