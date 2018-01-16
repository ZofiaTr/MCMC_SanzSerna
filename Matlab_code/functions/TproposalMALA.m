function prop = TproposalMALA(qNew, qOld,  dV, t)


     Tx = qNew - qOld + 0.5 * t *dV(qOld);
     
     Tabs = (Tx'*Tx ) / ( t);
     prop =   exp(-Tabs / 2.0);
end
