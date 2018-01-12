function [newrefseq,refIns_loc] = updateRef(dat,reffile,star_n)
    N = size(dat,1);
    part_I =cell(1,N); 
    for i=1:N
        part_I{i} = getInsertionPosition(dat(i));
    end
    part_I_all=cat(1,part_I{:});
    part_I_pos=part_I_all(any(part_I_all,2),:); %delete [0,0]
    [u,~,~]=unique(part_I_pos,'rows');
    refIns_loc=u;
    refIns_loc(:,1)=u(:,1)-star_n+1;% all the insertions
    if size(refIns_loc,1)==0
        newrefseq=reffile;
    else
        ins_piece = repmat('-',refIns_loc(1,2),1);
        newrefseq =[reffile(1:refIns_loc(1,1)-1), ins_piece.', reffile(refIns_loc(1,1):end)];

        for i=2:size(refIns_loc,1)
            ins_piece =repmat('-',refIns_loc(i,2),1);
            mid_loc = refIns_loc(i,1)+sum(refIns_loc(1:i-1,2));
            newrefseq = [newrefseq(1:mid_loc-1), ins_piece.', newrefseq(mid_loc:end)];
        end
    end
end


