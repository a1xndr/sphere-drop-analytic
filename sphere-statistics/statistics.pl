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
my $maxz=0.1;
while (my $line = <COORDS>) 
{
	chomp $line;
	my @words = split / /, $line;
	if($words[1]+$words[4]>$maxz){
	    $maxz = $words[4]+$words[1];
	}
	$volume+=$words[1]**3 * (4.0/3.0)*3.14159265;
	$count++;
}
print "\n    count:   ".$count;
print "\n    volume:  ".$volume;
print "\n    ratio:   ".$volume/(10*10*$maxz);
