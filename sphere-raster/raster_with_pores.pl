#!/usr/bin/perl -w
#
use warnings;
use GD::Simple;


@colors = {
"#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",

        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"
};

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
my $file2 = "pores";
open(COORDS2, $file2);
#binmode
while (my $line = <COORDS2>) 
{
	chomp $line;
	my @words = split / /, $line;
	push(@bubbles, { number=>$count, r => $words[0], x => $words[1], y => $words[2], z => $words[3], p=> $words[4]});
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
    $img->bgcolor('black');
    for(my $j=0; $j<=$count; $j++ ){
	if($spheres[$j]{z} + $spheres[$j]{r} >$layer && $spheres[$j]{z} - $spheres[$j]{r}<$layer){
	    $img->fgcolor('black');
	    $img->bgcolor('black');
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
    my $p;
    for(my $j=0; $j<=$bcount; $j++ ){
	if($bubbles[$j]{z} + $bubbles[$j]{r} >$layer && $bubbles[$j]{z} - $bubbles[$j]{r}<$layer){
	    $img->fgcolor('black');
            $p=$bubbles[$j]{p};
            $img->fgcolor($p%255, 255-$p%190, 255-$p%255);
            $img->bgcolor($p%255, 255-$p%190, 255-$p%255);
            #$img->bgcolor('white');
            #$img->fgcolor('white');
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
            #$img->string($bubbles[$j]{p});
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
