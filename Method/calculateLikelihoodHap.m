function [loglh] = calculateLikelihoodHap (read_alignseqi,quality_alignseqi,j,mapq,prim_scorek)
% alignmati: aligned reads and qscore at base i
% a: base type 1-5 'ATCG-'
% loglh:loglikelihood at base i
read_real = read_alignseqi(read_alignseqi~=-1);
qscore_real = quality_alignseqi(read_alignseqi~=-1);
mapq_real = mapq(read_alignseqi~=-1);
qscore_realnew = qscore_real;
qscore_realnew(read_real==5) = mapq_real(read_real==5);
qscore_j_rare = qscore_realnew(read_real == j);% store qscores of real reads that match j
qscore_nj_rare = qscore_realnew(read_real ~= j); % store qscores of real reads that do not match j
qscore_j = qscore_j_rare(all(qscore_j_rare,2),:);
qscore_nj = qscore_nj_rare(all(qscore_nj_rare,2),:);

if size(qscore_j,1) == 0
    loglh = NaN;
else  
    logseqerrprob_j = log10(1-10 .^ (-qscore_j / 10)); %log of seq error prob
    logseqerrprob_nj = -qscore_nj / 10+log10(0.25);
    loglh = sum(logseqerrprob_j)+sum(logseqerrprob_nj)+log(prim_scorek); %loglikelihood; prim_scorek is prior 
end
end


    
