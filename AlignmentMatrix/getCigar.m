function [sam_readi] =cigar(dati,sam_readsx,refIns_loc)

    %if S exist, remove S and corresponding informaiton from cigar information
    %and read
    cig_index_S = regexp(dati.CigarString,'S');

    cig_index_H = regexp(dati.CigarString,'H');

    if ~isempty(cig_index_H)
        str = dati.CigarString;
        expression = '\d*H';
        startIndex = regexp(str,expression,'match');
        dati.CigarString=regexprep(str,startIndex,'');
    end


    if ~isempty(cig_index_S)
        str = dati.CigarString;
        expression = '\d*S';
        startIndex = regexp(str,expression,'match');
        remainstr=str;
        sam_reads=sam_readsx;
        for t=1:size(startIndex,2)
            remainstr=strrep(remainstr, startIndex(t), '');
            cig_num_S =str2double(regexprep(startIndex{t},'[A-Z]',' '));
            if max(regexp(dati.CigarString,'M'))>cig_index_S(t)
                sam_reads=sam_reads(cig_num_S+1:end);
            else
                sam_reads=sam_reads(1:end-cig_num_S);
            end
        end
        newstr=remainstr{1};

    else
        newstr=dati.CigarString;
        sam_reads=sam_readsx;
    end

    cig_num =str2double(strsplit(regexprep(newstr,'[A-Z]',' ')));
    %find the index of match bases corresponding to reference
    cig_index = regexp(newstr,'[A-Z]');
    cig_index_M = regexp(newstr,'M');
    M_loc=zeros(1,length(cig_index_M));
    for t=1:length(cig_index_M)
        M_loc(t)=find(cig_index==cig_index_M(t));
    end
    loc_sami = dati.Position:dati.Position + sum(cig_num(M_loc))-1;

    %exclude deletion

    cig_index_D = regexp(newstr,'D');
    cig_num_exD=cig_num;
    cig_index_exD=cig_index;
    exD_loc=zeros(1,length(cig_index_D));
    for t=1:length(cig_index_D)
        exD_loc(t) =find(cig_index==cig_index_D(t));
    end
    cig_num_exD(exD_loc)=[];
    cig_index_exD(exD_loc)=[];

    %include index of insertion bases corresponding to reference

    cig_index_I = regexp(newstr,'I'); 
    I_real=zeros(length(cig_index_I),2);
    I_count=zeros(length(cig_index_I),1);
    for t=1:length(cig_index_I)
        sumI=sum(I_count(1:t));
        I_loc =find(cig_index_exD==cig_index_I(t));
        I_count(t) =cig_num_exD(I_loc); % number of base in the tth insertion
        I_index=loc_sami(sum(cig_num_exD(1:I_loc-1)));
        I_add = repmat(I_index+1,1,I_count(t));

        %find the location and length of insertion
        I_loc_real =find(cig_index==cig_index_I(t));
        I_index_real=dati.Position+sum(cig_num(1:I_loc_real-1))-sumI;
        I_real(t,:)=[I_index_real I_count(t)];%

        loc_sami=[loc_sami(1:sum(cig_num_exD(1:I_loc-1))) I_add loc_sami(sum(cig_num_exD(1:I_loc-1))+1:end)];
    end

    %include index of deletion bases corresponding to reference


    sam_readi=vertcat(loc_sami,sam_reads);

    for t=1:length(cig_index_D)
        D_loc =find(cig_index==cig_index_D(t));
        D_count =cig_num(D_loc); % number of base in the tth insertion
        D_index=sam_readi(1,sum(cig_num(1:D_loc-1)))+1;
        D_add=zeros(2,D_count);
        D_add(1,:)=D_index:D_index+D_count-1;
        D_add(2,:) = repmat(5,1,D_count);
        sam_readi(1,sum(cig_num(1:D_loc-1))+1:end)=sam_readi(1,sum(cig_num(1:D_loc-1))+1:end)+D_count;
        sam_readi=horzcat(sam_readi(:,1:sum(cig_num(1:D_loc-1))), D_add, sam_readi(:,sum(cig_num(1:D_loc-1))+1:end));
    end

    %including insertion from other cigarinformation
    totIns_loc = refIns_loc;
    for t = 1:size(I_real,1)
        tf=find(ismember(totIns_loc,I_real(t,:),'rows'));
        totIns_loc(tf,:)=[]; %delete insertions already in readi
    end

    for t=1:size(totIns_loc,1)
        Ins_loc=find(sam_readi(1,:)==totIns_loc(t,1),1);
        if Ins_loc>1
            Ins_add=zeros(2,totIns_loc(t,2));
            Ins_add(1,:) = repmat(sam_readi(1,Ins_loc),1,totIns_loc(t,2));
            Ins_add(2,:) = repmat(5,1,totIns_loc(t,2));
            sam_readi=[sam_readi(:,1:Ins_loc-1) Ins_add sam_readi(:,Ins_loc:end)];
        end
    end
end




