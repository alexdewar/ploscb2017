function acts=rkernmult(rkerns,im)
    acts = shiftdim(sum(sum(bsxfun(@times,rkerns,im),1),2),2);
end