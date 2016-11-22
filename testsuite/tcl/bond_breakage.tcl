
# Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
#  
# this file is part of espresso.
#  
# espresso is free software: you can redistribute it and/or modify
# it under the terms of the gnu general public license as published by
# the free software foundation, either version 3 of the license, or
# (at your option) any later version.
#  
# espresso is distributed in the hope that it will be useful,
# but without any warranty; without even the implied warranty of
# merchantability or fitness for a particular purpose.  see the
# gnu general public license for more details.
#  
# you should have received a copy of the gnu general public license
# along with this program.  if not, see <http://www.gnu.org/licenses/>. 

 

# This tests the bond breakage mechanism
# * UI
# * Bond breakage with tabulated bonded ia
source "tests_common.tcl"

require_feature "TABULATED"


# setup
setmd box_l 10 10 10

thermostat off
setmd time_step 0.01
setmd skin 0


# Interface, default state
if { "[bond_breakage]" != "off" } {
  error_exit "Bond breakage mechanism should be off by default" 
}

# Test adding handlers
bond_breakage add print_queue_entry break_simple_pair_bond
if { "[bond_breakage]" != "print_queue_entry break_simple_pair_bond" } {
  error_exit "Bond breakage mechanism should be off by default" 
}


# Test turning it off again
bond_breakage off
if { "[bond_breakage]" != "off" } {
  error_exit "Failed to disable the bond breakage mechanism" 
}


# Test bond breakage of tabulated bonds
# We are (mis)-using a tabulated lj as bonded tabulated potential

bond_breakage add print_queue_entry break_simple_pair_bond


proc test_for_bond { breakable_bond non_breakable_bond } {
part delete
eval "inter 0 $breakable_bond"
puts "inter 0 $breakable_bond"
part 0 pos 0 0 0
part 1 pos 1 0 0 bond 0 0
integrate 0 recalc_forces
puts [part 1 print bond]
if { "[part 1 print bond]" != "{ {0 0} } " } {
  error_exit "The bond was deleted, when it should not have been."
}
# Place the particle beyond the range of the potential
part 1 pos 1.5 0 0
integrate 0 recalc_forces
puts [part 1 print bond]
if { "[part 1 print bond]" != "{ } " } {
  error_exit "The bond is still there, even though it should have been removed"
}


# Now test that a broken bond error is thrown, when the tabulated bond is
# configured not to use the bond_breakage mechanism
eval "inter 0 $non_breakable_bond" 
puts [inter]
part 1 bond 0 0
if {![catch {integrate 0 recalc_forces} err]} {
    error_exit "no exception was thrown at over-stretched bond, although requested"
}
}


test_for_bond "harmonic 1 0 1.3 1" "harmonic 1 0 1.3"
test_for_bond "tabulated bond lj1.tab 1" "tabulated bond lj1.tab 0"

