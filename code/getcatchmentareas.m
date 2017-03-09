function ca=getcatchmentareas(views,ref_view_i,xg,yg)
    rmsdiffs = getRMSdiff(views,views(:,:,ref_view_i));
    heads = getIDFheads(xg,yg,shiftdim(rmsdiffs));
    ca = [getca(heads,ref_view_i,pi/4),getca(heads,ref_view_i,pi/2)];
end

function ca=getca(heads,ref_view_i,cutoff)
    success = heads<cutoff;
    success(ref_view_i) = true;
    bw = bwlabeln(success);
    ca = sum(bw(:)==bw(ref_view_i));
end