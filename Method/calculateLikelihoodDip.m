function [loglh] = calculateLikelihoodDip (read_alignseqi,quality_alignseqi,jv,mapq,prim_scorek)
% alignmati: aligned reads and qscore at base i
% a: base type 1-5 'ATCG-'
% loglh:loglikelihood at base i
read_real = read_alignseqi(read_alignseqi~=-1);
qscore_real = quality_alignseqi(read_alignseqi~=-1);
mapq_real = mapq(read_alignseqi~=-1);
qscore_realnew = qscore_real;
if isempty(qscore_real(read_real~=5))
qscore_realnew(read_real==5) = mapq_real(read_real==5)/4;
else
qscore_realnew(read_real==5) = mean(qscore_real(read_real~=5))
end
sublh=zeros(1,length(read_real));

if size(read_real(read_real==jv(1)),1)==0 && size(read_real(read_real==jv(2)),1)==0
    loglh=NaN;
else
    for subi=1:length(read_real)
        if read_real(subi)==jv(1) && read_real(subi)==jv(2)
            sublh(subi)=log10(1-10^(-qscore_realnew(subi)/ 10));
        elseif read_real(subi)~=jv(1) && read_real(subi)~=jv(2)
            sublh(subi) = -qscore_realnew(subi) / 10+log10(0.25);

        else    
            sublh(subi)=log10(0.5*(1-10^(-qscore_realnew(subi)/10))+0.125*10^(-qscore_realnew(subi)/10));
        end
    end
    loglh=sum(sublh)+log(prim_scorek);

end
readsiz=size(read_real,1);
if sum(jv==[5,5])==2 && sum(read_real==5)>readsiz*0.8
     loglh=0;
end
end

    
    

