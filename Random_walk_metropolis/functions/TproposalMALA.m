function prop = TproposalMALA(pNew, pOld,  dV, t)


     Tx = pNew - pOld + t.*dV(pOld);
   
     Tabs = (Tx'*Tx ) / (2 * t);
     prop =   exp(-Tabs / 2.0);
end
