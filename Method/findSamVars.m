function [varsamfile,varsamfilesup]=findSamVars(hap_base, base_state,refseq,aln_refindex,diploid,method)
k = 1;
supk=1;
if diploid==1
r = sum(base_state~=1 & base_state~=-1 & base_state~=10);
varsam=ones(r,3);
 if method==1 
    for i=1:length(hap_base)
       if base_state(i) >= 2 && base_state(i) <= 5
           varsam(k,:)=[aln_refindex(i) refseq(aln_refindex(i)) hap_base(i)];
           k=k+1;
       elseif base_state(i)> 20
           varsam(k,:)=[aln_refindex(i) 5 hap_base(i)];
           k=k+1;
       end 
    end 
 elseif method==2
    for i=1:length(hap_base)
       if base_state(i) >= 2 && base_state(i) <= 5
           varsam(k,:)=[aln_refindex(i) refseq(aln_refindex(i)) hap_base(i)];
           k=k+1;
       elseif base_state(i)> 5 && base_state(i) < 10
           varsam(k,:)=[aln_refindex(i) 5 hap_base(i)];
           k=k+1;
       end 
    end 
 end
  
    ref_var=convertSequence(varsam(:,2), 'ACGT-');
    hap_var=convertSequence(varsam(:,3), 'ACGT-');
    if r==1 % the special case when there is only one snp or indel
         ref_var=[ref_var 'A'];
         hap_var=[hap_var 'A'];
         varsamfile=dataset([varsam(:,1);0],ref_var.',hap_var.');
         varsamfile(double(varsamfile(:,1))==0,:)=[];
    else
         varsamfile=dataset(varsam(:,1),ref_var.',hap_var.');
    end
    varsamfilesup=[];
   
elseif diploid==2
         
     if method==1
         r = sum(base_state~=1 & base_state~=-1);
         varsam=zeros(r,4);
         for i=1:size(hap_base,2)
             if base_state(i) >= 2 && base_state(i) <= 15
                  varsam(k,:)=[aln_refindex(i) refseq(aln_refindex(i)) hap_base(1,i) hap_base(2,i)];
                  k=k+1;
             elseif base_state(i)> 20
                  varsam(k,:)=[aln_refindex(i) 5 hap_base(1,i) hap_base(2,i)];
                  k=k+1;
             end 
         end 
     elseif method==2 
         r = sum(base_state~=1 & base_state~=-1 & base_state~=30);
         varsam=zeros(r,4);
         for i=1:size(hap_base,2)
             if base_state(i) >= 2 && base_state(i) <= 15
                varsam(k,:)=[aln_refindex(i) refseq(aln_refindex(i)) hap_base(1,i) hap_base(2,i)];
                k=k+1;
             elseif base_state(i)> 15 && base_state(i) < 30
                varsam(k,:)=[aln_refindex(i) 5 hap_base(1,i) hap_base(2,i)];
                k=k+1;
             end 
         end  
     end
   
     clvarsam=zeros(r,3);
     clvarsam(:,1:2)=varsam(:,1:2);
     
     for i=1:r
        if varsam(i,3)==varsam(i,4)
            clvarsam(i,3)=varsam(i,3);
        else
           if varsam(i,3)==varsam(i,2)
               clvarsam(i,3)=varsam(i,4);
           elseif varsam(i,4)==varsam(i,2)
               clvarsam(i,3)=varsam(i,3);
           else
               supclvarsam(supk,:)=[varsam(i,1) varsam(i,2) varsam(i,3) varsam(i,4)];
               supk=supk+1;
           end
        end
     end
     
      clvarsam(sum((clvarsam==0),2)>0,:) = [];
      ref_var=convertSequence(clvarsam(:,2), 'ACGT-');
      hap_var1=convertSequence(clvarsam(:,3), 'ACGT-');
      
     if size(clvarsam,1)==1
         ref_varsup=[ref_var 'A'];
         hap_var1sup=[hap_var1 'A'];
         varsamfile=dataset([clvarsam(:,1);0],ref_varsup.',hap_var1sup.');
         varsamfile(double(varsamfile(:,1))==0,:)=[];
     else
         varsamfile=dataset(clvarsam(:,1),ref_var.',hap_var1.');
     end


         if supk>1
             ref_varsup=convertSequence(supclvarsam(:,2), 'ACGT-');
             hap_var1sup=convertSequence(supclvarsam(:,3), 'ACGT-');
             hap_var2sup=convertSequence(supclvarsam(:,4), 'ACGT-');
             if supk==2
                  ref_varsup=[ref_varsup 'A'];
                  hap_var1sup=[hap_var1sup 'A'];
                  hap_var2sup=[hap_var2sup 'A'];
                  varsamfilesup=dataset([supclvarsam(:,1);0],ref_varsup.',hap_var1sup.',hap_var2sup.');
                  varsamfilesup(double(varsamfilesup(:,1))==0,:)=[];
             else
                  varsamfilesup=dataset(supclvarsam(:,1),ref_varsup.',hap_var1sup.',hap_var2sup.');
             end  
         else
             varsamfilesup=[]; 
         end    
end
end

    
    
    



