function [dat,refseq_tru,aln_read,aln_quality,aln_refindex]=getAlignMatrixMain(samfilename,start_n,reffile)


file.mainDir = pwd;

file.samDir = 'data';

file.samFile = fullfile(file.mainDir, file.samDir,samfilename);

dat = samread(file.samFile); 
%transform to A/C/G/T/N to 1/2/3/4/7
refseq_tru = [];
refseq_tru(reffile=='A') = 1; 
refseq_tru(reffile=='C') = 2;
refseq_tru(reffile=='G') = 3;
refseq_tru(reffile=='T') = 4;
refseq_tru(reffile=='N') = 7;
[newrefseq,refIns_loc] = updateRef(dat,reffile,start_n);% renwew ref based on sam file
[aln_read,aln_quality] = updateSam(dat,newrefseq,refIns_loc,start_n); %align read to newrefseq based on sam file

%index ref
ref_align=aln_read(1,:);
aln_refindex = cumsum(ref_align ~= 5);

end
