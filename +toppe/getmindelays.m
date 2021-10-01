function d = getmindelays(system)

d.rf.pre = system.start_core_rf + system.myrfdel;
d.daq.pre = system.start_core_daq + system.daqdel;
d.grad.pre = system.start_core_grad;
d.post = system.timetrwait + system.timessi; 
