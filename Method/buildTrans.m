function [transrow] = buildTrans(tprobi,hetrate)
transrow = zeros(1,30);
transrow(30) = tprobi(4) * hetrate/32;
transrow(1) = tprobi(1) * (1-hetrate)*(1-transrow(30));
transrow(2:4) = tprobi(2) * (1-hetrate)/3*(1-transrow(30));
transrow(5:7) = (tprobi(1) + tprobi(2)/3) * hetrate/4*(1-transrow(30));
transrow(8) = (tprobi(1) + tprobi(3)) * hetrate/4*(1-transrow(30));
transrow([9:10,12]) = tprobi(2) * hetrate/6*(1-transrow(30));
transrow([11,13,14]) = (tprobi(2)/3 + tprobi(3)) * hetrate/4*(1-transrow(30));
transrow(15) = tprobi(3) * (1-hetrate)*(1-transrow(30));
transrow(16:19) = tprobi(4) * (1-hetrate)/4*(1-transrow(30));
transrow([20:22,24,25,27]) = tprobi(4) * hetrate/8*(1-transrow(30));
transrow([23,26,28,29]) = tprobi(4) * hetrate/16*(1-transrow(30));

end


