function d = getmindelays(system)

d.rf.pre = system.toppe.start_core_rf + system.toppe.myrfdel;
d.daq.pre = system.toppe.start_core_daq + system.toppe.daqdel;
d.grad.pre = system.toppe.start_core_grad;
d.post = system.toppe.timetrwait + system.toppe.timessi; 
