function [hap_base] = findDipBaseStat(base_state, ref_state)
n=length(base_state);
hap_base=zeros(2,n)-1;
hap_base(1,base_state==1)=ref_state(1,base_state==1);
hap_base(2,:)=hap_base(1,:);
state_vecA = [1 1;2 2;3 3;4 4;1 2;1 3;1 4;1 5;2 3;2 4;2 5;3 4;3 5;4 5;5 5];
state_vecC = [2 2;1 1;3 3;4 4;1 2;2 3;2 4;2 5;1 3;1 4;1 5;3 4;3 5;4 5;5 5];
state_vecG = [3 3;1 1;2 2;4 4;1 3;2 3;3 4;3 5;1 2;1 4;1 5;2 4;2 5;4 5;5 5];
state_vecT = [4 4;1 1;2 2;3 3;1 4;2 4;3 4;4 5;1 2;1 3;1 5;2 3;2 5;3 5;5 5];
state_vec = [1 1;2 2;3 3;4 4;1 2;1 3;1 4;1 5;2 3;2 4;2 5;3 4;3 5;4 5;5 5];

for i=1:n

  if base_state(i) > 1
    if ref_state(i)==1
        hap_temp = vertcat(state_vecA, state_vec);
    elseif ref_state(i)==2
        hap_temp = vertcat(state_vecC, state_vec);
    elseif ref_state(i)==3
        hap_temp = vertcat(state_vecG, state_vec);
    elseif ref_state(i)==4
        hap_temp = vertcat(state_vecT, state_vec);
    else
        hap_temp = vertcat(state_vecA, state_vec);  
    end
    hap_base(:,i)=hap_temp (base_state(i),:).';
  end
end
end

