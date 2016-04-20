#!/usr/bin/perl -w
#author: Chen Zhaoqiang
#contact: 744891290@qq.com


use strict;
use Data::Dumper;

require "rotate.pl";
print $ARGV[0],"\n";
my $jsontext;
{
local $/="";
open FH,$ARGV[0];
$jsontext=<FH>;
$jsontext=~s/(?<="):/=>/g;            #"
}

my $parahash=eval($jsontext);


print $parahash->{'xyz_file'},"\n";
print join "-",@{$parahash->{'rotate_ids'}},"\n";

my @rotatebond=@{$parahash->{'rotate_bond'}};

my %coorhash;

open FH,$parahash->{'xyz_file'} or die "can't find the file $parahash->{'xyz_file'} ";
my $atom_count=0;
while(<FH>)
{
	if ($_=~/\S+/)
	{
	
		my @elements=split(/\s+/,$_);
		$atom_count++;
		print "$atom_count  AAA\n";
		$coorhash{$atom_count}{'atomtype'}=$elements[0];
		$coorhash{$atom_count}{'x'}=$elements[1];
		$coorhash{$atom_count}{'y'}=$elements[2];
		$coorhash{$atom_count}{'z'}=$elements[3];
		
    }
}

#print Dumper(\%coorhash);

my @graph=@{$parahash->{'rotate_ids'}}; #'

my $diheangle;
$diheangle=90/180 * 3.14159;#现在二面角的角度 修改

my $gjfheaders='# PM6 scan nosymm

S2

0 1
';
my %outpu;
#&rotatee(\%coorhash,$diheangle,\@graph,\@rotatebond);

%outpu=%{dclone(&rotatee(\%coorhash,$diheangle,\@graph,\@rotatebond))};

for(1..$atom_count)
{
	#print "$outpu{$_}{'atomtype'}   $outpu{$_}{'x'}     $outpu{$_}{'y'}   $outpu{$_}{'z'}\n";
	my $x=sprintf("%.8f",$outpu{$_}{'x'});
	my $y=sprintf("%.8f",$outpu{$_}{'y'});
	my $z=sprintf("%.8f",$outpu{$_}{'z'});
	printf("$outpu{$_}{'atomtype'} %+15.8s  %+15.8s %+15.8s\n",$x,$y,$z);
	
}

my @rotaterange=@{$parahash->{'rotate_range'}};
print join "---",@rotaterange;

#die ("test");
my $dirname; 
($dirname=$parahash->{'xyz_file'})=~s/\.xyz/_output/;
mkdir($dirname);
#die "$rotaterange[1]";
for (my $diheangle=$rotaterange[0];$diheangle<=$rotaterange[1];$diheangle+=$rotaterange[2])
{
	#die ($diheangle);
	my  $filename=$dirname.'_'.$diheangle.'.gjf';
	
	#die($filename);
	print $filename,"BBB\n";
    my $dihe2=$diheangle/180 * 3.14159;#现在二面角的角度 修改

	my %outpu;
#&rotatee(\%coorhash,$diheangle,\@graph,\@rotatebond);

    %outpu=%{dclone(&rotatee(\%coorhash,$dihe2,\@graph,\@rotatebond))};
   
    
	my $outgjffile='./'."$dirname".'/'.$filename;
	my $fh;
	open $fh,">$outgjffile";
    &writegjf(\%outpu,$atom_count,$fh);
	
}


sub writegjf
{
	my %outpu=%{$_[0]};
	my $atom_count=$_[1];
	my $fh=$_[2];
	print $fh $gjfheaders;
	for(1..$atom_count)
	{
		#print "$outpu{$_}{'atomtype'}   $outpu{$_}{'x'}     $outpu{$_}{'y'}   $outpu{$_}{'z'}\n";
		my $x=sprintf("%.8f",$outpu{$_}{'x'});
		my $y=sprintf("%.8f",$outpu{$_}{'y'});
		my $z=sprintf("%.8f",$outpu{$_}{'z'});
		my $line=sprintf("$outpu{$_}{'atomtype'} %+15.8s  %+15.8s %+15.8s\n",$x,$y,$z);
		print $fh $line;
	}
	
	
	
}
