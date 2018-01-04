function prop = TproposalMALA(pNew, pOld,  dV, t)


     Tx = pNew - pOld + 0.5 * t *dV(pOld);
     
     Tabs = (Tx'*Tx ) / ( t);
     prop =   exp(-Tabs / 2.0);
end
