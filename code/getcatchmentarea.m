function ca=getcatchmentarea(views,ref_view_i,xg,yg)
    rmsdiffs = getRMSdiff(views,views(:,:,ref_view_i));
    heads = getIDFheads(xg,yg,shiftdim(rmsdiffs));
    success = heads<pi/2;
    success(ref_view_i) = true;
    bw = bwlabeln(success);
    ca = sum(bw(:)==bw(ref_view_i));
end