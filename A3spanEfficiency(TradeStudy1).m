function spanEff = spanEfficiency(AR,sweepLE)
% (Raymer 12.49)
spanEff = 4.61.*(1 - 0.045.*AR.^0.68).*cosd(sweepLE).^0.15 - 3.1;
end