function in = isdentro(x,y,idx_bordo_mant, idx_bordo_osta,vertices)
    [xm,ym] = poly2ccw(vertices(1,idx_bordo_mant),vertices(2,idx_bordo_mant));
    [xo,yo] = poly2cw(vertices(1,idx_bordo_osta),vertices(2,idx_bordo_osta));
    
    xv = [xm, NaN, xo];
    yv = [ym, NaN, yo];
    in = inpolygon(x,y,xv,yv);
end