#!/usr/bin/perl -w

# Rounding procedure of the mighty PerlMonks...
# Halleluja!

sub round {
    my( $num, $prec )= @_;
    return int( $num/$prec + 0.5 - ($num<0) ) * $prec;
}

# Procedure to calculate the distance between 2 Vectors
sub abstand {
    my( $a1,$a2,$a3,$b1,$b2,$b3)= @_;
    return ( (($a1-$b1)**2 + ($a2-$b2)**2 + ($a3-$b3)**2)**0.5);
}

# Check input
if (@ARGV<7) {
    die "\n\n
###########################################################################
                               INFLATEGRO 
                  Written by Christian Kandt, (c) 2005-2007

   Kandt C, Ash WL, Tieleman DP (2007): Setting up and running molecular
        dyanmics simulations of membrane proteins. Methods 41:475-488
###########################################################################

USAGE: 
-----
INFLATEGRO  bilayer.gro  scaling_factor  lipid_residue_name  cutoff   inflated_bilayer.gro   gridsize areaperlipid.dat (protein)
\n\n";
}

if (!(open (INPUT, $ARGV[0]))) {
    print "Eeeeek! No $ARGV[0] at all!\n\n";
    die;
}

# Read lipid coordinates
$scale    = $ARGV[1];
$name     = $ARGV[2];
$cutoff   = $ARGV[3]*0.1;
$output   = $ARGV[4];
$gridsize = $ARGV[5];
$area     = $ARGV[6];

$switch = 0;
if (@ARGV==8 and $ARGV[7] eq "protein") {
    print "Well, just the protein then....\n";
    $switch=1;
}

$zaehler      = 1;
$counter      = 1;
$protein_xmax = -1000;
$protein_ymax = -1000;
$protein_xmin = 1000;
$protein_ymin = 1000;

print STDOUT "Reading..... \n";
my $line_count = 0;

while (<INPUT>) {
    chomp;
    my $line = $_;
    $line_count++;

    # Saltar las dos primeras líneas (.gro header y número de átomos)
    next if ($line_count <= 2);

    # Detectar la última línea (dimensiones de la caja)
    if ($line =~ /^\s*([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s*$/) {
        $box_x = $1;
        $box_y = $2;
        $box_z = $3;
        next;
    }

    # Procesar átomos por columnas fijas según el estándar GROMACS
    if (length($line) >= 44) {
        my $resnum_str  = substr($line, 0, 5);
        my $resname_str = substr($line, 5, 5);
        my $atmname_str = substr($line, 10, 5);
        my $atmnum_str  = substr($line, 15, 5);
        my $x_str       = substr($line, 20, 8);
        my $y_str       = substr($line, 28, 8);
        my $z_str       = substr($line, 36, 8);

        # LIMPIEZA TOTAL: Quitar espacios al inicio y al final de cada campo extraído
        $resnum_str  =~ s/^\s+|\s+$//g;
        $resname_str =~ s/^\s+|\s+$//g;
        $atmname_str =~ s/^\s+|\s+$//g;
        $atmnum_str  =~ s/^\s+|\s+$//g;
        $x_str       =~ s/^\s+|\s+$//g;
        $y_str       =~ s/^\s+|\s+$//g;
        $z_str       =~ s/^\s+|\s+$//g;

        # Clasificación
        if ($resname_str eq $name) {
            $resnum_l[$zaehler]  = $resnum_str;
            $resname_l[$zaehler] = $resname_str;
            $atmname_l[$zaehler] = $atmname_str;
            $atmnum_l[$zaehler]  = $atmnum_str;
            $x_l[$zaehler]       = $x_str;
            $y_l[$zaehler]       = $y_str;
            $z_l[$zaehler]       = $z_str;
            $zaehler++;
        } elsif ($resname_str ne "SOL") {
            $resnum_p[$counter]  = $resnum_str;
            $resname_p[$counter] = $resname_str;
            $atmname_p[$counter] = $atmname_str;
            $atmnum_p[$counter]  = $atmnum_str;
            $x_p[$counter]       = $x_str;
            $y_p[$counter]       = $y_str;
            $z_p[$counter]       = $z_str;

            if ($x_p[$counter] > $protein_xmax) { $protein_xmax = $x_p[$counter]; }
            if ($x_p[$counter] < $protein_xmin) { $protein_xmin = $x_p[$counter]; }
            if ($y_p[$counter] > $protein_ymax) { $protein_ymax = $y_p[$counter]; }
            if ($y_p[$counter] < $protein_ymin) { $protein_ymin = $y_p[$counter]; }
            $counter++;
        }
    }
}
close (INPUT);

$zaehler--;
$counter--;
$totalatmn = $zaehler + $counter;

# Converting nm into A
$protein_xmin = $protein_xmin*10;
$protein_xmax = $protein_xmax*10;
$protein_ymin = $protein_ymin*10;
$protein_ymax = $protein_ymax*10;

# New boxsize
my $old_box_x = $box_x;
my $old_box_y = $box_y;
$box_x = $box_x * $scale;
$box_y = $box_y * $scale;

# 1. Determinar el centro de masa (COM) de la proteína PRIMERO
$xcenter = 0.5 * $box_x;
$ycenter = 0.5 * $box_y;
$xcom = 0; 
$ycom = 0;
$translatex_p = 0;
$translatey_p = 0;

if ($counter == 0) {
    print "No protein coordinates found. Using center of old box...\n";
    $xcom = 0.5 * $old_box_x;
    $ycom = 0.5 * $old_box_y;
}
if ($counter > 0) {
    print STDOUT "Centering protein....\n";
    $xpsum = 0; $ypsum = 0;
    for ($k=1; $k<=$counter; $k++) {
         $xpsum = $xpsum + $x_p[$k];
         $ypsum = $ypsum + $y_p[$k];
    }
    $xcom  = $xpsum / $counter;
    $ycom  = $ypsum / $counter;
    $translatex_p = $xcenter - $xcom;
    $translatey_p = $ycenter - $ycom;
}

# 2. ESCALADO DE LÍPIDOS Y CONTEO DE PCOUNT (Ahora va ANTES de la determinación de hojuelas)
print STDOUT "Scaling lipids radially....\n";
$pcount = 1;

for ($k=1; $k<=$zaehler; $k++) {
    if ($atmname_l[$k] =~ /^P/) {
        $pxneu = $xcenter + ($x_l[$k] - $xcom) * $scale;
        $pyneu = $ycenter + ($y_l[$k] - $ycom) * $scale;
        
        $res = $resnum_l[$k];
        $translatex_l[$res] = $pxneu - $x_l[$k];
        $translatey_l[$res] = $pyneu - $y_l[$k];
        $phosz[$pcount] = $z_l[$k];
        $pcount++;
    }
}

$pcount--;

if ($pcount == 0) {
    die "\n¡ERROR FATAL!: El script no encontró ningún átomo que empiece con 'P' en tu residuo lipídico '$name'.\n".
        "Por favor, revisa cómo se llama el átomo de fósforo en tu archivo .gro (ej. P, P8, PO4).\n\n";
}

print "There are $pcount lipids...\n";
$atomperlipid = $zaehler/$pcount;
print "with $atomperlipid atoms per lipid..\n";

# 3. DETERMINACIÓN DE HOJUELAS SUPERIOR E INFERIOR
print "\nDetermining upper and lower leaflet...\n";
$middle = 0;
for ($p=1; $p<=$pcount; $p++) {
    $middle = $middle + $phosz[$p];
}
$middle = $middle / $pcount;

$uppercount = 0;
$lowercount = 0;
$upper = 0;
$lower = 0;

for ($p=1; $p<=$pcount; $p++) {
    if ($phosz[$p] > $middle) {
        $upper = $upper + $phosz[$p];
        $uppercount++;
    }
    if ($phosz[$p] < $middle) {
        $lower = $lower + $phosz[$p];
        $lowercount++;
    }
}
$upper = $upper / $uppercount;
$lower = $lower / $lowercount;

print "$uppercount lipids in the upper...\n";
print "$lowercount lipids in the lower leaflet \n\n"; 

# Checking for protein lipid overlap
$upper_rm = 0;
$lower_rm = 0;

if ($cutoff > 0 && $switch == 0) {
    print "Checking for overlap....\n";
    print "...this might actually take a while....\n";

    $overlapcount = 0;
    for ($k=1; $k<=$zaehler; $k++) {
        $uppercheck = 0;
        $lowercheck = 0;
        $progress = ($k/$zaehler)*100;
        $progress = round ($progress,2);
        print STDOUT "$progress % done...\r";

        if ($atmname_l[$k] =~ /^P/) {
            $res = $resnum_l[$k];
            $overlap[$res] = 0;
             
            for ($i=1; $i<=$counter; $i++) {
                if ($atmname_p[$i] eq "CA") {
                    $distance = abstand($x_l[$k]+$translatex_l[$res], $y_l[$k]+$translatey_l[$res], $z_l[$k], $x_p[$i]+$translatex_p, $y_p[$i]+$translatey_p, $z_p[$i]);
                    if ($distance <= $cutoff) {
                        $overlap[$res] = 1;
                        if ($z_l[$k] > $middle) { $uppercheck = 1; }
                        if ($z_l[$k] < $middle) { $lowercheck = 1; }
                    }
                }
            }
            $overlapcount = $overlapcount + $overlap[$res];
            if ($uppercheck == 1) { $upper_rm++; }
            if ($lowercheck == 1) { $lower_rm++; }
        }
    }
    print "\nThere are $overlapcount lipids within cut-off range...\n";
    print "$upper_rm will be removed from the upper leaflet...\n";
    print "$lower_rm will be removed from the lower leaflet...\n\n";
}

$newlipids     = $pcount - $upper_rm - $lower_rm;
$newupper      = $uppercount - $upper_rm;
$newlower      = $lowercount - $lower_rm;
$totalatmn_new = $totalatmn - ($upper_rm + $lower_rm)*$atomperlipid;

# Writing scaled bilayer & centered protein
print STDOUT "Writing scaled bilayer & centered protein...\n";
open(OUTPUT, ">$output") or die "Cannot write to $output";

print OUTPUT "InflateGRO reordered and fixed. Box scaled by $scale\n";
print OUTPUT "$totalatmn_new\n";

for ($k=1; $k<=$counter; $k++) {
    $newx = $x_p[$k] + $translatex_p;
    $newy = $y_p[$k] + $translatey_p;
    printf OUTPUT "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", $resnum_p[$k], $resname_p[$k], $atmname_p[$k], $atmnum_p[$k], $newx, $newy, $z_p[$k]; 
}

if ($switch == 0) {
    for ($k=1; $k<=$zaehler; $k++) {
        $res  = $resnum_l[$k];
        $newx = $x_l[$k] + $translatex_l[$res];
        $newy = $y_l[$k] + $translatey_l[$res];

        if ($cutoff > 0) {
            if (defined $overlap[$res] && $overlap[$res] == 0) {
                printf OUTPUT "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", $resnum_l[$k], $resname_l[$k], $atmname_l[$k], $atmnum_l[$k], $newx, $newy, $z_l[$k]; 
            }
        } else {
            printf OUTPUT "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", $resnum_l[$k], $resname_l[$k], $atmname_l[$k], $atmnum_l[$k], $newx, $newy, $z_l[$k]; 
        }
    }
}

printf OUTPUT "%10.5f%10.5f%10.5f\n", $box_x, $box_y, $box_z;
close OUTPUT;

# Calculate Area per lipid
print "\n\nCalculating Area per lipid...\n";

$protein_xmax = int ($protein_xmax+1);
$protein_xmin = int ($protein_xmin);
$protein_ymax = int ($protein_ymax+1);
$protein_ymin = int ($protein_ymin);

$xrange = $protein_xmax - $protein_xmin;
$yrange = $protein_ymax - $protein_ymin;

print "Protein X-min/max: $protein_xmin    $protein_xmax\n";
print "Protein Y-min/max: $protein_ymin    $protein_ymax\n";
print "X-range: $xrange A    Y-range: $yrange A\n";

if ($protein_xmin != 0 or $protein_xmax != 0){
    for ($k=1; $k<=$counter; $k++) {
        $x_p[$k] = 10 * $x_p[$k] - $protein_xmin;
        $y_p[$k] = 10 * $y_p[$k] - $protein_ymin;
    }   
}

# Building 2D grid on protein coordinates
print "Building $xrange X $yrange 2D grid on protein coordinates...\n";
for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
        $grid[$x][$y] = 0;
    }
}

# Calculating area occupied by protein
print "Calculating area occupied by protein..\n";
print "full TMD..\n";

for ($k=1; $k<=$counter; $k++) {
    if ($z_p[$k] >= $lower and $z_p[$k] <= $upper) {
        $x = int( $x_p[$k] / $gridsize);  
        $y = int( $y_p[$k] / $gridsize);
        $grid[$x][$y] = 1 if ($x >= 0 && $y >= 0 && $x <= $xrange/$gridsize && $y <= $yrange/$gridsize);
    }
    $progress = $k/$counter * 100;
    $progress = round ($progress,1);
    print "$progress % done...\r";
}

$howmany = 0;
for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
         $howmany = $howmany + $grid[$x][$y];
    }
}

$areaprotein_total = ($gridsize)**2 * $howmany * 0.01;
$arealipid_total   = ($box_x * $box_y - $areaprotein_total) / ($newlipids * 0.5);

print "upper TMD..\n";
for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
        $grid[$x][$y] = 0;
    }
}

for ($k=1; $k<=$counter; $k++) {
    if ($z_p[$k] >= $middle and $z_p[$k] <= $upper) {
        $x = int( $x_p[$k] / $gridsize);  
        $y = int( $y_p[$k] / $gridsize);
        $grid[$x][$y] = 1 if ($x >= 0 && $y >= 0 && $x <= $xrange/$gridsize && $y <= $yrange/$gridsize);
    }
    $progress = $k/$counter * 100;
    $progress = round ($progress,1);
    print "$progress % done...\r";
}

$howmany = 0;
for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
         $howmany = $howmany + $grid[$x][$y];
    }
}

$areaprotein_upper = ($gridsize)**2 * $howmany * 0.01;
$arealipid_upper   = ($box_x * $box_y - $areaprotein_upper) / ($newupper);

print "lower TMD..\n";
for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
        $grid[$x][$y] = 0;
    }
}

for ($k=1; $k<=$counter; $k++) {
    if ($z_p[$k] >= $lower and $z_p[$k] <= $middle) {
        $x = int( $x_p[$k] / $gridsize);  
        $y = int( $y_p[$k] / $gridsize);
        $grid[$x][$y] = 1 if ($x >= 0 && $y >= 0 && $x <= $xrange/$gridsize && $y <= $yrange/$gridsize);
    }
    $progress = $k/$counter * 100;
    $progress = round ($progress,1);
    print "$progress % done...\r";
}

$howmany = 0;
for ($x=0; $x<=$xrange/$gridsize; $x=$x+1) {
    for ($y=0; $y<=$yrange/$gridsize; $y=$y+1) {
         $howmany = $howmany + $grid[$x][$y];
    }
}

$areaprotein_lower = ($gridsize)**2 * $howmany * 0.01;
$arealipid_lower   = ($box_x * $box_y - $areaprotein_lower) / ($newlower);

print "Area per protein: $areaprotein_total nm^2\n";
print "Area per lipid: $arealipid_total nm^2\n\n";
print "Area per protein, upper half: $areaprotein_upper nm^2\n";
print "Area per lipid, upper leaflet : $arealipid_upper nm^2\n\n";
print "Area per protein, lower half: $areaprotein_lower nm^2\n";
print "Area per lipid, lower leaflet : $arealipid_lower nm^2\n\n";

print STDOUT "Writing Area per lipid...\n";
open(AREAOUT, ">$area") or die "Cannot write to $area";
print AREAOUT "$arealipid_total \n";
close AREAOUT;

print "Done!\n\n\n";
