proc save_vtf {fname} {
 set f [open $fname w]
 writevsf $f
 writevcf $f
 close $f
}

proc save_pov {filename} {
  global nx
  global nz
  set f [open $filename w]
  for {set i 0} {$i <[setmd n_part]} {incr i} {
#  for {set i 0} {$i <[expr ($nx+$nz)+1]} {incr i} {
#    if{[part $i print type]>0}{
#      continue}
     puts $f [part $i print pos quat]}
  }
#}

#originaly:
setmd skin 0.4
#setmd min_global_cut 1.5 

setmd box_l 50 50 50
setmd time_step 0.01
thermostat langevin 1 1 

set nx 10 
set nz 0 

set aspect 1.5
set sigma1 1
set sigma2 1
set sigma3 1.5 
set epsilon1 1
set epsilon2 1
set epsilon3 1 
set steps 10

# Location of the lj minimum 
set lj_fac [expr pow(2,1./6.)]

#setmd min_global_cut [expr $sigma3]
setmd min_global_cut 1.5 

inter 0 0 gay-berne $epsilon1 $epsilon2 $epsilon3 $sigma1 $sigma2 $sigma3 2.5 2 1 
inter 0 harmonic 0 1.5
inter 1 angle_harmonic 800 [PI]
inter 3 harmonic 800 0

part 0 pos 50 50 50 

for {set i 1} {$i <$nx} {incr i} {
  part [setmd n_part] pos [expr $sigma3 * $i+1] 0 0 quat 1 0 0 0 mass 1
}

for {set i 1} {$i <$nz} {incr i} {
  part [setmd n_part] pos 0 0 [expr $sigma3 * $i+1 ] quat 1 0 0 0 mass 1 
}

puts "The observable gives the distance to the specified particle"
#puts "[analyze distto 0]"
puts "[analyze mindist]"
#on_collision bind_at_point_of_collision [expr $sigma3 *$aspect *$lj_fac *1.05] 0 1 1
on_collision triangle_binding [expr $sigma3 *$lj_fac *1.05] 0 1 1 0.15
save_vtf "before_integrate.vtf"
integrate 0
save_vtf "after_integrate0.vtf" 
#puts("Total number of particles after integrate 0 is:")
puts [setmd n_part]
puts "Observable analyzes distance to the specified particle"
puts "[analyze distto 0]"

for {set i 1} {$i <=$steps} {incr i} {
  integrate 1
  puts [setmd n_part]
  if {fmod($i,100)==0} {
    save_pov "bond_info-triaxial_ellipsoids_t_$i.dat"
  } 
  if {fmod($i,100)==0} {
    save_vtf "trixial_t_$i.vtf"
  }
}

save_pov "bond_info-triaxial_ellipsoids.dat"
save_vtf "after_integrate_steps.vtf" 
puts [setmd n_part]

exit
