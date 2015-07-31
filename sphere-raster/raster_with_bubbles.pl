#!/usr/bin/perl -w
#
use warnings;
use GD::Simple;


my $file = $ARGV[0] or die

#create a new image (width, height)
my $img = GD::Simple->new(200, 100);

my $count =0;

my @AoH;
#local $/ = undef;
#open FILE, $ARGV[$1] or die "Couldn't open file: $!";
open(COORDS, $file);
my @captures;
#binmode
my $volume=0; 
while (my $line = <COORDS>) 
{
	chomp $line;
	my @words = split / /, $line;
	push(@spheres, { number=>$count, r => $words[0], x => $words[1], y => $words[2], z => $words[3] });
	$volume+=$words[0]**3 * (4/3)*3.14159265;
	$count++;
}

my $bcount = 0;
my $file2 = "bubbles";
open(COORDS2, $file2);
#binmode
while (my $line = <COORDS2>) 
{
	chomp $line;
	my @words = split / /, $line;
	push(@bubbles, { number=>$count, r => $words[0], x => $words[1], y => $words[2], z => $words[3] });
	$volume+=$words[0]**3 * (4/3)*3.14159265;
	$bcount++;
}

use Data::Dumper qw(Dumper);;   
print Dumper \@spheres;

my $MIN=0;
my $MAX=10;
my $STEP=0.03;
my $layer = $MIN;

my $i=0;
my $name = "slice";
while($layer<10 ){
    my $num;
    my $img = GD::Simple->new(1000, 1100);
    my $area = 0;

    $img->font('Arial');
    $img->fontsize(10);
    $img->bgcolor('black');
    $img->fgcolor('black');
    $img->rectangle(0,0,1000,1020);
    $img->fgcolor('white');
    $img->bgcolor('white');
    for(my $j=0; $j<=$count; $j++ ){
	if($spheres[$j]{z} + $spheres[$j]{r} >$layer && $spheres[$j]{z} - $spheres[$j]{r}<$layer){
	    $img->fgcolor('white');
	    $img->bgcolor('white');
	    $img->moveTo(abs(1000*($spheres[$j]{x}))/$MAX,abs(1000*(10-$spheres[$j]{y}))/$MAX);
	    #print 500*($spheres[$j]{x})/$MAX."\n";
            #print $j.'\n';
	    $diameter=sqrt(($spheres[$j]{r})**2-($spheres[$j]{z}-$layer)**2)*2;
	    $area += 3.14159265*($diameter/2)**2;
	    $diameter=(1000/$MAX)*$diameter;
	    $img->ellipse($diameter,$diameter);
	    $img->moveTo(abs(1000*($spheres[$j]{x}))/$MAX,abs(1000*(10-$spheres[$j]{y}))/$MAX);
	    $img->fgcolor('red');
	    $img->bgcolor('red');
            #$img->string($j);
            #print "r:".$spheres[$j]{r}." x:".$spheres[$j]{x}." y:".$spheres[$j]{y}." z:".$spheres[$j]{z} ."";
	    $num++;
	}
    }
    for(my $j=0; $j<=$bcount; $j++ ){
	if($bubbles[$j]{z} + $bubbles[$j]{r} >$layer && $bubbles[$j]{z} - $bubbles[$j]{r}<$layer){
	    $img->fgcolor('black');
	    $img->bgcolor('red');
	    $img->moveTo(abs(1000*($bubbles[$j]{x}))/$MAX,abs(1000*(10-$bubbles[$j]{y}))/$MAX);
	    #print 500*($bubbles[$j]{x})/$MAX."\n";
            print $j.'\n';
	    $diameter=sqrt(($bubbles[$j]{r})**2-($bubbles[$j]{z}-$layer)**2)*2;
	    $area += 3.14159265*($diameter/2)**2;
	    $diameter=(1000/$MAX)*$diameter;
	    $img->ellipse($diameter,$diameter);
	    $img->moveTo(abs(1000*($bubbles[$j]{x}))/$MAX,abs(1000*(10-$bubbles[$j]{y}))/$MAX);
	    $img->fgcolor('white');
	    $img->bgcolor('white');
            #$img->string($j);
            print "r:".$bubbles[$j]{r}." x:".$bubbles[$j]{x}." y:".$bubbles[$j]{y}." z:".$bubbles[$j]{z} ."";
	    $num++;
	}
    }
    my $area_ratio=$area/100;
    $area_ratio = sprintf("%.2f", $area_ratio);
    $area = sprintf("%.2f", $area);
    $img->bgcolor('white');
    $img->fgcolor('white');
    $img->rectangle(0,1000,1000,1050);
    $img->bgcolor('black');
    $img->fgcolor('black');
    $img->moveTo(10,1045);
    $img->fontsize(37);
    $img->font('Arial:bold');
    $level = sprintf("%.2f", $layer);
    $img->string('z: '.$level);
    $img->moveTo(500,1045);
    $img->string('num visible: '.$num);
    $img->moveTo(10,1095);
    $img->string('Sphere Area: '.$area);
    $img->moveTo(500,1095);
    $img->string('Area Fraction: '.$area_ratio);
    $img->fontsize(20);
    open my $out, '>', 'img'.$i.'.png' or die;
    binmode $out;
    print $out $img->png;
    $i++;
    $layer= $layer+$STEP;
}

print "    volume:".$volume;
print "    ratio:".$volume/1000;
