function [varsam,varsamsup] = runVarCalling(aln_read,aln_quality,dat,aln_refindex,refseq_tru,diploid,prim,tprob,hetrate)

 nread=size(dat,1);
 mapq=zeros(1,nread);
 for i = 1 : nread
     mapq(i)= dat(i).MappingQuality;
 end

read_align=aln_read;
read_align(1,:)=[];
read_alignseq=read_align;
quality_alignseq=aln_quality;
ref_align=aln_read(1,:);

if diploid==1
    [base_state] = findStateHap(read_alignseq, quality_alignseq,mapq,ref_align,prim,tprob);
    hap_base = findHapBaseStat(base_state, aln_read(1,:));
elseif  diploid==2
    [base_state] = findStateDip(read_alignseq, quality_alignseq,mapq,ref_align,prim,tprob,hetrate);
    hap_base = findDipBaseStat(base_state, aln_read(1,:));
end

[varsam,varsamsup] = findSamVars(hap_base, base_state,refseq_tru,aln_refindex,diploid,2);
end


