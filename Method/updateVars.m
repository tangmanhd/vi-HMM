function [varout_sort ]= updateVars(varmat,varmatsup,reffile,start_n)

%%find all SNPs 
SNP_loc=find( double(varmat(:,2))~=45 & double(varmat(:,3))~=45);
if size(SNP_loc,1)==0
    VarNames = {'Var1','Var2','Var3'};
    varmat_SNP=dataset([],[],[],'VarNames',VarNames);
else
    varmat_SNP=varmat(SNP_loc,:);
end


%%find all Inserations;
ins_loc=find( double(varmat(:,2))==45 & double(varmat(:,3))~=45);
if size(ins_loc,1)==0
    varmat_ins=[];
else
    varmat_ins_n=varmat(ins_loc,:);
    ins_chr_uniq=double(unique(varmat_ins_n(:,1)));
    varmat_ins=cell(size(ins_chr_uniq,1),3);
    for i = 1:size(ins_chr_uniq,1)
        vari=varmat_ins_n(double(varmat_ins_n(:,1))==ins_chr_uniq(i),:);
        var2i=reffile(ins_chr_uniq(i));
        cellvari=dataset2cell(vari(:,3));
        var3i=strjoin(cellvari(2:end),''); %first cell is var3, the colname.thus start from 2:end
        varmat_ins{i,1}=ins_chr_uniq(i);
        varmat_ins{i,2}=var2i;
        varmat_ins{i,3}=strcat(var2i,var3i);   
    end
end

%%find all delation;
del_loc=find(double(varmat(:,2))~=45 & double(varmat(:,3))==45);
if size(del_loc,1)==0
    varmat_del=[];
else
    varmat_del_n=varmat(del_loc,:);
    del_chr_uniq=double(unique(varmat_del_n(:,1)));
    varmat_del=cell(size(del_chr_uniq,1),3);
    for i = 1:size(del_chr_uniq,1)
        vari=varmat_del_n(double(varmat_del_n(:,1))==del_chr_uniq(i),:);
        cellvari=dataset2cell(vari(:,2));
        var2i=strjoin(cellvari(2:end),''); %first cell is var3, the colname.thus start from 2:end
        var3i=reffile(del_chr_uniq(i)-1);
        varmat_del{i,1}=del_chr_uniq(i)-1;
        varmat_del{i,2}=strcat(var3i,var2i);
        varmat_del{i,3}=var3i;   
    end
end

%%%for supplementary var:base has different alleles
if size(varmatsup,2)==0 %%the case that no varmatsup file
    varmatsup_SNP=[];
    varmatsup_SNP_del=[];
    varmatsup_ins=[];
else 
    SNPsup_loc=find(double(varmatsup(:,2))~=45 & double(varmatsup(:,3))~=45 & double(varmatsup(:,4))~=45);
    if size(SNPsup_loc,1)==0
        varmatsup_SNP=[];
    else
        varmatsup_SNP_n=varmatsup(SNPsup_loc,:);
        varmatsup_SNP=cell(size(SNPsup_loc,1),3);
        for i = 1:size(SNPsup_loc,1)
            varmatsup_SNP{i,1}=double(varmatsup_SNP_n(i,1));
            var2i=dataset2cell(varmatsup_SNP_n(i,2));
            var3i=dataset2cell(varmatsup_SNP_n(i,3));
            var4i=dataset2cell(varmatsup_SNP_n(i,4));
            varmatsup_SNP{i,2}=var2i(2);
            varmatsup_SNP{i,3}=strcat(var3i(2),{','},var4i(2));   
        end
    end

    SNP_delsup_loc=find( double(varmatsup(:,2))~=45 & double(varmatsup(:,3))~=45 & double(varmatsup(:,4))==45);
    if size(SNP_delsup_loc,1)==0
        varmatsup_SNP_del=[];
    else
        varmatsup_SNP_del_n=varmatsup(SNP_delsup_loc,:);
        varmatsup_SNP_del=cell(size(SNP_delsup_loc,1),3);
        for i = 1:size(SNP_delsup_loc,1)
            varmatsup_SNP_del{i,1}=double(varmatsup_SNP_del_n(i,1))-1;
            var3i=dataset2cell(varmatsup_SNP_del_n(i,3));
            varmatsup_SNP_del{i,2}=reffile(varmatsup_SNP_del{i,1}:varmatsup_SNP_del{i,1}+1);
            varmatsup_SNP_del{i,3}=strcat(strcat(reffile(varmatsup_SNP_del{i,1}),var3i(2)),{','},reffile(varmatsup_SNP_del{i,1}));   
        end
    end

    inssup_loc=find( double(varmatsup(:,2))==45 & double(varmatsup(:,3))~=45 & double(varmatsup(:,4))~=45);
    if size(inssup_loc,1)==0
        varmatsup_ins=[];
    else
        varmatsup_ins_n=varmatsup(inssup_loc,:);
        inssup_chr_uniq=double(unique(varmatsup_ins_n(:,1)));
        varmatsup_ins=cell(size(inssup_chr_uniq,1),3);
        for i = 1:size(inssup_chr_uniq,1)
            vari=varmatsup_ins_n(double(varmatsup_ins_n(:,1))==inssup_chr_uniq(i),:);
            var2i=reffile(inssup_chr_uniq(i));
            cellvari=dataset2cell(vari(:,3));
            var3i=strjoin(cellvari(2:end),'');
            cellvar2i=dataset2cell(vari(:,4));
            var4i=strjoin(cellvar2i(2:end),'');
            varmatsup_ins{i,1}=inssup_chr_uniq(i);
            varmatsup_ins{i,2}=var2i;
            varmatsup_ins{i,3}=strcat(strcat(var2i,var3i),{','},strcat(var2i,var4i));  
        end
    end
end
varout=[dataset2cell(varmat_SNP);varmat_ins; varmat_del; varmatsup_SNP;varmatsup_SNP_del;varmatsup_ins];
ds = cell2dataset(varout);
org_loc=double(ds(:,1))+start_n-1;
org_loc_table=mat2dataset(org_loc,'VarNames','org_loc');
ds_new=horzcat(org_loc_table,ds);
varout_sort = sortrows(ds_new(:,[1 3 4]));
