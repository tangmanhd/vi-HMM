function [emissionmatrix] = buildEmissionMatrixHap(ref_sec,read_sec,quality_sec,mapq,prim)
logemissionmatrix=zeros(10,size(ref_sec,2));
if size(prim,1)==0
    prim=zeros(2,5)+1;
end
for k=1:size(ref_sec,2)
    a = 1:5;
    M = ref_sec(1,k);
    S1 = min(a(a ~= M));
    S2 = min(a(a ~= M & a ~= S1));
    S3 = min(a(a ~= M & a ~= S1 & a ~= S2));
    S4 = max(a(a ~= M));

    hap_temp = [M,S1,S2,S3,S4,1,2,3,4,5];
    if ref_sec(1,k)==5 %check insert
        logemissionmatrix(1:5,k) = NaN;
        for s=6:10
            logemissionmatrix(s,k) = calculateLikelihoodHap(read_sec(:,k),quality_sec(:,k),hap_temp(s),mapq,prim(2,s-5));
        end
    else
        for s=1:5
            logemissionmatrix(s,k) = calculateLikelihoodHap(read_sec(:,k),quality_sec(:,k),hap_temp(s),mapq,prim(1,s));
        end
        logemissionmatrix(6:10,k) = NaN;
    end   
end
emissionmatrix=logemissionmatrix;
end

