#!/usr/bin/perl -w

#author: Chen Zhaoqiang
#contact: 744891290@qq.com



use strict;
use Storable qw(dclone);
use Data::Dumper;
	



unless(caller)
{

###input START#####
my %coor;
my @graph;
my $diheangle;
my @rotatebond=(3,8); #8是要旋转的平面!!! 修改
while(<DATA>)
{
	
	my @elements=split(/\s+/,$_);
	
	$coor{$.}{'atomtype'}=$elements[0];
	$coor{$.}{'x'}=$elements[1];
	$coor{$.}{'y'}=$elements[2];
	$coor{$.}{'z'}=$elements[3];
	
	
}
@graph=(8,34,35,36);
$diheangle=179.344/180 * 3.14159;#现在二面角的角度 修改

#print $diheangle,"\n";
###input  END#####


###output  START#####
my %outpu;
%outpu=%{dclone(&rotatee(\%coor,$diheangle,\@graph,\@rotatebond))};

#print Dumper(\%outpu);	
for(1..47)
{
	print "$outpu{$_}{'atomtype'}   $outpu{$_}{'x'}     $outpu{$_}{'y'}   $outpu{$_}{'z'}\n";
	
}
######output END#####

}
#旋转轴的方向指向要旋转的平面，旋转的角度是二面角的角度取反
sub rotatee
{
	my %hash=%{dclone($_[0])};
	my $rotateangle=-$_[1];
	my @graph=@{$_[2]};
	my @rotatebond=@{$_[3]};
	print join "bond",@rotatebond,"\n";
    print join "->",@graph,"\n";
  #  print Dumper(\%hash);
	my @A=($hash{$rotatebond[0]}{'x'},$hash{$rotatebond[0]}{'y'},$hash{$rotatebond[0]}{'z'});
	my @B=($hash{$rotatebond[1]}{'x'},$hash{$rotatebond[1]}{'y'},$hash{$rotatebond[1]}{'z'});#旋转的平面
	print join "|",@B,"\n";
	my $dist=($B[0]-$A[0])**2+($B[1]-$A[1])**2+($B[2]-$A[2])**2;
	$dist=sqrt($dist);
	
	print "lenght of bond is $dist\n";
	my @direction=(($B[0]-$A[0])/$dist,($B[1]-$A[1])/$dist,($B[2]-$A[2])/$dist);


	my $T=
	[
		[1,0,0,0],
		[0,1,0,0],
		[0,0,1,0],
		[-$A[0],-$A[1],-$A[2],1]
	];

	my $TF=
	[
		[1,0,0,0],
		[0,1,0,0],
		[0,0,1,0],
		[$A[0],$A[1],$A[2],1]


	];
	

	my @B2=@direction;

	my $cosU=$B2[2]/sqrt($B2[1]**2+$B2[2]**2);
	my $sinU=$B2[1]/sqrt($B2[1]**2+$B2[2]**2);

	my $RX=
	[
		[1,0,0,0],
		[0,$cosU,$sinU,0],
		[0,-$sinU,$cosU,0],
		[0,0,0,1]

	];
	
	my $RXF=
	[
		[1,0,0,0],
		[0,$cosU,-$sinU,0],
		[0,$sinU,$cosU,0],
		[0,0,0,1]
	];


	my $cosV=sqrt($B2[1]**2+$B2[2]**2)/sqrt($B2[0]**2+$B2[1]**2+$B2[2]**2);
	my $sinV=-$B2[0]/sqrt($B2[0]**2+$B2[1]**2+$B2[2]**2);

	my $RY=
	[

		[$cosV,0,-$sinV,0],
		[0,1,0,0],
		[$sinV,0,$cosV,0],
		[0,0,0,1]
		
	];
		
	my $RYF=
	[
		[$cosV,0,$sinV,0],
		[0,1,0,0],
		[-$sinV,0,$cosV,0],
		[0,0,0,1]

	];

	my $cosW=cos($rotateangle);
	my $sinW=sin($rotateangle);

	my $RZ=
	[
		[$cosW,$sinW,0,0],
		[-$sinW,$cosW,0,0],
		[0,0,1,0],
		[0,0,0,1]


	];

  foreach my $id(@graph)
	{
		my @X=($hash{$id}{'x'},$hash{$id}{'y'},$hash{$id}{'z'});
		push @X,1;
		my $vect=
		[
			[@X]
		];
		
		
		my $result=&mmmult($vect,$T,$RX,$RY,$RZ,$RYF,$RXF,$TF);
    #my $result=pdl($vect) x pdl($T) x pdl($RX) x pdl($RY) x pdl($RZ) x pdl($RYF) x pdl($RXF) x pdl($TF);
		 #print "$result\n";
		my $resultvect=join "|",map(@{$_},@{$result});
		


		if($resultvect=~/(\S+)\|(\S+)\|(\S+)\|/)
		{
			$hash{$id}{'x'}=$1;
			$hash{$id}{'y'}=$2;
			$hash{$id}{'z'}=$3;


		}



	}

	return \%hash;


}









sub mmmult
{
	my @parameters=@_;
  
	
  my $result=shift @parameters;

	foreach my $para(@parameters)
	{
		

		$result=&mmult($result,$para);
		
		
	}
	
	return $result;

}
		




#################function4 matrix mult


sub mmult {

    my ($m1,$m2) = @_;

    my ($m1rows,$m1cols) = matdim($m1);

    my ($m2rows,$m2cols) = matdim($m2);


    unless ($m1cols == $m2rows) {  # raise exception


        die "IndexError: matrices don't match: $m1cols != $m2rows";


    }


    my $result = [  ];


    my ($i, $j, $k);


    for $i (range($m1rows)) {


        for $j (range($m2cols)) {


            for $k (range($m1cols)) {


                $result->[$i][$j] += $m1->[$i][$k] * $m2->[$k][$j];


            }


        }


    }


    return $result;


}


sub range { 0 .. ($_[0] - 1) }


sub veclen {


    my $ary_ref = $_[0];
  


    my $type = ref $ary_ref;


    if ($type ne "ARRAY") { die "$type is bad array ref for $ary_ref" }


    return scalar(@$ary_ref);


}



sub matdim {


    my $matrix = $_[0];


    my $rows = veclen($matrix);


    my $cols = veclen($matrix->[0]);


    return ($rows, $cols);


}



__DATA__
N                2.372  -0.007   0.043  0.00  0.00             
N               -0.775   4.726  -0.592  0.00  0.00             
C                3.952  -2.460  -0.386  0.00  0.00             
C                2.593  -2.466   0.310  0.00  0.00             
C                4.594  -1.085  -0.218  0.00  0.00             
C                1.713  -1.332  -0.199  0.00  0.00             
C                3.664   0.014  -0.718  0.00  0.00             
C                4.851  -3.547   0.174  0.00  0.00             
C                1.528   1.163  -0.126  0.00  0.00             
C                0.237   1.262   0.429  0.00  0.00             
C                2.074   2.333  -0.691  0.00  0.00             
C               -0.504   2.437   0.289  0.00  0.00             
C                1.350   3.505  -0.807  0.00  0.00             
C                0.030   3.561  -0.340  0.00  0.00             
C               -0.424   0.234   1.299  0.00  0.00             
C               -1.599  -0.489   0.720  0.00  0.00             
C               -1.745  -0.607  -0.665  0.00  0.00             
C               -2.557  -1.054   1.566  0.00  0.00             
C               -1.823   5.032   0.397  0.00  0.00             
C               -2.835  -1.279  -1.202  0.00  0.00             
C               -3.652  -1.725   1.035  0.00  0.00             
C               -3.785  -1.834  -0.347  0.00  0.00             
O               -0.066   0.058   2.449  0.00  0.00             
Cl              -5.109  -2.647  -0.994  0.00  0.00           
H                3.797  -2.651  -1.477  0.00  0.00             
H                2.728  -2.388   1.408  0.00  0.00             
H                2.082  -3.433   0.135  0.00  0.00             
H                5.549  -1.039  -0.778  0.00  0.00             
H                4.858  -0.916   0.845  0.00  0.00             
H                1.492  -1.468  -1.282  0.00  0.00             
H                0.734  -1.378   0.315  0.00  0.00             
H                4.152   1.003  -0.576  0.00  0.00             
H                3.490  -0.093  -1.811  0.00  0.00             
H                5.038  -3.408   1.248  0.00  0.00             
H                4.401  -4.540   0.045  0.00  0.00             
H                5.826  -3.556  -0.330  0.00  0.00             
H                3.113   2.338  -1.058  0.00  0.00             
H               -1.526   2.490   0.702  0.00  0.00             
H                1.821   4.383  -1.265  0.00  0.00             
H               -0.201   5.523  -0.764  0.00  0.00             
H               -0.991  -0.162  -1.328  0.00  0.00             
H               -2.442  -0.973   2.655  0.00  0.00             
H               -1.443   5.487   1.324  0.00  0.00             
H               -2.361   4.109   0.671  0.00  0.00             
H               -2.538   5.727  -0.061  0.00  0.00             
H               -2.952  -1.373  -2.287  0.00  0.00             
H               -4.406  -2.169   1.694  0.00  0.00             
