function [hap_base]=findHapBaseStat(base_state, ref_state)
n=length(base_state);
hap_base=zeros(1,n)-1;
hap_base(base_state==1)=ref_state(1,base_state==1);
 a = 1:5;
for i=1:n
    if base_state(i)==2
        newa=a(a~=ref_state(1,i));
        hap_base(i)=newa(1);
    elseif base_state(i)==3
        newa=a(a~=ref_state(1,i));
        hap_base(i)=newa(2);
    elseif base_state(i)==4
        newa=a(a~=ref_state(1,i));
        hap_base(i)=newa(3);
    elseif base_state(i)==5
        newa=a(a~=ref_state(1,i));
        hap_base(i)=newa(4);
    elseif base_state(i)==6  
        hap_base(i)=1; 
    elseif base_state(i)==7  
        hap_base(i)=2;
    elseif base_state(i)==8  
        hap_base(i)=3;
    elseif base_state(i)==9  
        hap_base(i)=4;  
    elseif base_state(i)==10
        hap_base(i)=-1; 
    end
end
end


        
