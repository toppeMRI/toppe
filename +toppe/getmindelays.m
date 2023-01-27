function d = getmindelays(system)

d.rf.pre = system.start_core_rf + system.psd_rf_wait;
d.daq.pre = system.start_core_daq + system.psd_grd_wait;
d.grad.pre = system.start_core_grad;
d.post = system.timetrwait + system.timessi; 
