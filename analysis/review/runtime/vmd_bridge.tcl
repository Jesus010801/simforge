# SimForge ↔ VMD bridge — protocol simforge-vmd/1
#
# Static script.  Every dynamic value (port, token, file paths, formats)
# arrives through environment variables and is only ever used as a Tcl
# *variable value* — never spliced into code, never evaluated.  Received
# lines are parsed by an allow-list; there is no eval / subst / uplevel.
#
# VMD here only displays coordinates and reports its frame.  It computes
# nothing, transforms nothing and writes no trajectory.
#
# Environment:
#   SIMFORGE_VMD_PORT, SIMFORGE_VMD_TOKEN
#   SIMFORGE_VMD_STRUCTURE, SIMFORGE_VMD_STRUCTURE_TYPE        (gro | pdb)
#   SIMFORGE_VMD_TRAJECTORY, SIMFORGE_VMD_TRAJECTORY_TYPE      (xtc | trr)
#
# Frames: the structure file's own coordinates are deleted after loading, so
# VMD frame i == display frame i of the ReviewDataset (both zero-based).

namespace eval ::simforge {
    variable protocol "simforge-vmd/1"
    variable maxline 512
    variable sock ""
    variable state "init"
    variable molid -1
    variable nframes 0
    variable pending_seq ""
    variable pending_frame ""
}

# ── pure helpers (tested with plain tclsh) ───────────────────────────────────

proc ::simforge::is_uint {text} {
    return [regexp {^(0|[1-9][0-9]{0,9})$} $text]
}

# One SimForge → VMD line  →  {WELCOME} | {GOTO seq frame} | {BYE} | {ERROR code}
proc ::simforge::parse_command {line} {
    variable protocol
    variable maxline
    if {[string length $line] > $maxline} { return [list ERROR line_too_long] }
    if {![regexp {^[ -~]*$} $line]} { return [list ERROR bad_line] }
    set parts [split $line " "]
    set cmd [lindex $parts 0]
    switch -exact -- $cmd {
        WELCOME {
            if {[llength $parts] == 2 && [lindex $parts 1] eq $protocol} { return [list WELCOME] }
            return [list ERROR bad_welcome]
        }
        GOTO {
            if {[llength $parts] != 3} { return [list ERROR bad_goto] }
            set seq [lindex $parts 1]
            set fr [lindex $parts 2]
            if {![is_uint $seq] || ![is_uint $fr]} { return [list ERROR bad_goto] }
            return [list GOTO $seq $fr]
        }
        BYE {
            if {[llength $parts] == 1} { return [list BYE] }
            return [list ERROR bad_bye]
        }
        default { return [list ERROR unknown_command] }
    }
}

# Report a frame; it acknowledges the pending command iff it is that frame.
proc ::simforge::frame_report {frame} {
    variable pending_seq
    variable pending_frame
    if {$pending_seq ne "" && $frame == $pending_frame} {
        set ack $pending_seq
        set pending_seq ""
        set pending_frame ""
    } else {
        set ack NONE
    }
    return "FRAME $frame $ack"
}

proc ::simforge::clean {text} {
    regsub -all {[^ -~]} $text " " text
    return [string range $text 0 200]
}

# ── VMD / socket side ────────────────────────────────────────────────────────

proc ::simforge::send {line} {
    variable sock
    if {$sock eq ""} { return }
    if {[catch {puts $sock $line; flush $sock}]} { ::simforge::shutdown }
}

proc ::simforge::shutdown {} {
    variable sock
    variable state
    set state "closed"
    if {$sock ne ""} {
        catch {fileevent $sock readable {}}
        catch {close $sock}
        set sock ""
    }
    after idle quit
}

proc ::simforge::on_frame {name1 name2 op} {
    variable molid
    variable state
    global vmd_frame
    if {$state ne "ready" || $name2 != $molid} { return }
    ::simforge::send [::simforge::frame_report $vmd_frame($molid)]
}

proc ::simforge::goto {seq frame} {
    variable nframes
    variable pending_seq
    variable pending_frame
    if {$frame >= $nframes} {
        ::simforge::send "ERROR frame_out_of_range $frame"
        return
    }
    set pending_seq $seq
    set pending_frame $frame
    animate goto $frame
    # no trace write (already on that frame): acknowledge explicitly, once
    if {$pending_seq ne ""} { ::simforge::send [::simforge::frame_report $frame] }
}

proc ::simforge::load {} {
    variable molid
    variable nframes
    variable state
    ::simforge::send LOADING
    set st $::env(SIMFORGE_VMD_STRUCTURE)
    set stype $::env(SIMFORGE_VMD_STRUCTURE_TYPE)
    set tr $::env(SIMFORGE_VMD_TRAJECTORY)
    set ttype $::env(SIMFORGE_VMD_TRAJECTORY_TYPE)
    if {[catch {
        set m [mol new $st type $stype waitfor all]
        set own [molinfo $m get numframes]
        mol addfile $tr type $ttype waitfor all molid $m
        if {$own > 0} { animate delete beg 0 end [expr {$own - 1}] $m }
    } err]} {
        ::simforge::send "ERROR load_failed [::simforge::clean $err]"
        return
    }
    set molid $m
    set nframes [molinfo $m get numframes]
    trace add variable ::vmd_frame($m) write ::simforge::on_frame
    set state "ready"
    ::simforge::send "READY $m [molinfo $m get numatoms] $nframes [vmdinfo version]"
}

proc ::simforge::handle {line} {
    variable state
    set msg [::simforge::parse_command $line]
    switch -exact -- [lindex $msg 0] {
        WELCOME {
            if {$state eq "hello"} {
                set state "loading"
                after idle ::simforge::load
            } else {
                ::simforge::send "ERROR unexpected_welcome"
            }
        }
        GOTO {
            if {$state ne "ready"} {
                ::simforge::send "ERROR not_ready"
            } else {
                ::simforge::goto [lindex $msg 1] [lindex $msg 2]
            }
        }
        BYE { ::simforge::shutdown }
        default { ::simforge::send "ERROR [lindex $msg 1]" }
    }
}

proc ::simforge::on_read {} {
    variable sock
    if {$sock eq ""} { return }
    if {[catch {gets $sock line} n]} { ::simforge::shutdown; return }
    if {$n < 0} {
        if {[eof $sock]} { ::simforge::shutdown; return }
        # an unterminated line may not grow without bound
        if {[chan pending input $sock] > 4096} { ::simforge::shutdown }
        return
    }
    ::simforge::handle $line
}

proc ::simforge::main {} {
    variable sock
    variable state
    variable protocol
    set sock [socket 127.0.0.1 $::env(SIMFORGE_VMD_PORT)]
    fconfigure $sock -buffering line -translation lf -encoding ascii -blocking 0
    fileevent $sock readable ::simforge::on_read
    set state "hello"
    ::simforge::send "HELLO $protocol $::env(SIMFORGE_VMD_TOKEN)"
}

if {![info exists ::simforge_test_mode]} { ::simforge::main }
