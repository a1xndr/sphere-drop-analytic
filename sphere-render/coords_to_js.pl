#!/usr/bin/perl -w
#
use warnings;

my $file = $ARGV[0] or die


my $count =0;

my @AoH;
#local $/ = undef;
#open FILE, $ARGV[$1] or die "Couldn't open file: $!";
open(COORDS, $file);
my @captures;
#binmode
print " var coords = [ ";
while (my $line = <COORDS>) 
{
        chomp $line;
	my @words = split / /, $line;
        if(@words==5){
        print "     ,[".$words[1].",".$words[2].",".$words[3].",".$words[4]."]";
	$count++;
    }
    }
print " ];";
