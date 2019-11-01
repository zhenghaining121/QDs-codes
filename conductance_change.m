function conductance = conductance_change(uplimit,downlimit,bins,cond)

conductance = uplimit-cond/bins*(uplimit-downlimit);


