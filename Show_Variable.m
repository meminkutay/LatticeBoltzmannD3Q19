function Show_Variable(Var, SI, SJ)
    Varmatrix = reshape(Var,[SI,SJ]);
    cla,imagesc(Varmatrix), daspect([1,1,1])
%     colorbar
    impixelinfo