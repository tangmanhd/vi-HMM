function [emissionmatrix] = buildEmissionMatrixDip(ref_sec,read_sec,quality_sec,mapq,prim)
logemissionmatrix=zeros(30,size(ref_sec,2));
if size(prim,1)==0
    prim=zeros(2,15)+1;
end
state_vecA = [1 1;2 2;3 3;4 4;1 2;1 3;1 4;1 5;2 3;2 4;2 5;3 4;3 5;4 5;5 5];
state_vecC = [2 2;1 1;3 3;4 4;1 2;2 3;2 4;2 5;1 3;1 4;1 5;3 4;3 5;4 5;5 5];
state_vecG = [3 3;1 1;2 2;4 4;1 3;2 3;3 4;3 5;1 2;1 4;1 5;2 4;2 5;4 5;5 5];
state_vecT = [4 4;1 1;2 2;3 3;1 4;2 4;3 4;4 5;1 2;1 3;1 5;2 3;2 5;3 5;5 5];

state_vec = [1 1;2 2;3 3;4 4;1 2;1 3;1 4;1 5;2 3;2 4;2 5;3 4;3 5;4 5;5 5];

for k=1:size(ref_sec,2)

    if ref_sec(1,k)==1
        hap_temp = vertcat(state_vecA, state_vec);
    elseif ref_sec(1,k)==2
        hap_temp = vertcat(state_vecC, state_vec);
    elseif ref_sec(1,k)==3
        hap_temp = vertcat(state_vecG, state_vec);
    elseif ref_sec(1,k)==4
        hap_temp = vertcat(state_vecT, state_vec);
    end
    
    if ref_sec(1,k)==5 %check insert
        logemissionmatrix(1:15,k) = NaN;
        for s=16:30
            logemissionmatrix(s,k) = calculateLikelihoodDip(read_sec(:,k),quality_sec(:,k),state_vec(s-15,:),mapq,prim(2,s-15));
        end
    else
        for s=1:15
            logemissionmatrix(s,k) = calculateLikelihoodDip(read_sec(:,k),quality_sec(:,k),hap_temp(s,:),mapq,prim(1,s));
        end
        logemissionmatrix(16:30,k) = NaN;
    end   
end
emissionmatrix=logemissionmatrix;
end

