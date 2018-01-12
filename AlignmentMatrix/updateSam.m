function [aln_read,aln_quality]= updateSam(dat,newrefseq,refIns_loc,start_n)

    % convert 'ATCG-' to '12345'
    sam_refseq(newrefseq=='A') = 1;
    sam_refseq(newrefseq=='C') = 2;
    sam_refseq(newrefseq=='G') = 3;
    sam_refseq(newrefseq=='T') = 4;
    sam_refseq(newrefseq=='-') = 5;
    sam_refseq(newrefseq=='N') = 7;
    sam_ix = cumsum(sam_refseq ~= 5);
    aln_read = zeros(size(dat,1),length(newrefseq))-1;  %matrix to store the aligned reads

    for i=1:size(dat,1)
        sam_read = dat(i).Sequence;
        sam_readsx=zeros(1,length(sam_read));
        %convert 'ATCG-' to '12345'
        sam_readsx(sam_read=='A') = 1;
        sam_readsx(sam_read=='C') = 2;
        sam_readsx(sam_read=='G') = 3;
        sam_readsx(sam_read=='T') = 4;
        sam_readsx(sam_read=='-') = 5;
        dat(i).Position = dat(i).Position - start_n + 1;
        sam_readi = getCigar(dat(i),sam_readsx,refIns_loc); %readi with all the cigar information and refIns_loc
        readi_start=find(sam_ix==sam_readi(1,1),1,'first'); %start location of readi in ref
        readi_end=readi_start+size(sam_readi,2)-1; %end location of readi in ref
        aln_read(i,readi_start:readi_end)=sam_readi(2,:);
    end
    aln_read=vertcat(sam_refseq,aln_read);

    aln_quality = zeros(size(dat,1),length(newrefseq))-1;

    for i=1:size(dat,1)
        sam_quality = dat(i).Quality;
        sam_qualityx = double(sam_quality)-33; %convert char to number
        sam_qualityi = getCigar(dat(i),sam_qualityx,refIns_loc);
        qualityi_start=find(sam_ix==sam_qualityi(1,1),1,'first');
        qualityi_end=qualityi_start+size(sam_qualityi,2)-1;
        aln_quality(i,qualityi_start:qualityi_end)=sam_qualityi(2,:);
    end
end











    


