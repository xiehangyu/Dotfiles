import m5
from m5.objects import *
system=System()
system.clk_domain=SrcClockDomain()
system.clk_domain.clock='1GHz'
system.clk_domain.voltage_domain=VoltageDomain()
system.mem_mode='timing'
system.mem_ranges=[AddrRange("512MB")]
system.cpu=X86TimingSimpleCPU()
system.membus=SystemXBar()
system.cpu.icache_port=system.membus.cpu_side_ports
system.cpu.dcache_port=system.membus.cpu_side_ports
system.cpu.createInterruptController()