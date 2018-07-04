set filename "./Draw_list.dat"
set listA {}
set listB {}
set listC {}
set f [open $filename r]
foreach line [split [read $f] \n] {
    lappend listA [lindex $line 0]
    lappend listB  [lindex $line 1]
    lappend listC  [lindex $line 2]
    }

draw delete all
draw material Glass1
if {[llength $listA] == [llength $listB]} {
   for {set x 0} {$x < [llength $listA]} {incr x} {
	set c [lindex $listC $x]
        draw color $c
        set residA [lindex $listA $x]
        set residB [lindex $listB $x]
#       echo $residA $residB
        set sel1 [atomselect top "resid $residA and name CA"]
        set sel2 [atomselect top "resid $residB and name CA"]
        lassign [$sel1 get {x y z}] start_d
        lassign [$sel2 get {x y z}] end_d
        draw cylinder $start_d $end_d radius 0.3 resolution 100 filled yes
       }
}


