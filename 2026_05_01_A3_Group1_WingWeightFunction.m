function W_w = wingWeight(Kdw,Kvs,Wdg,Nz,Sw,AR,t_c,taper,sweepLE,S_csw)
% (Raymer 15.1)
W_w = (0.0103.*Kdw.*Kvs.*((Wdg.*Nz).^0.5).*(Sw.*0.622)*(AR.^0.785) ...
    .*(t_c.^-0.4).*((1 + taper).^0.05).*(cosd(sweepLE).^-1).*(S_csw.^0.04));
end

