#! /usr/bin/perl -w

#
# Author: Alfredo Velasco
# The Bioinformatics Toolsmith, The University of Tulsa
#

# This utility module will download epigenetics files given a link and a file with
# the

#package Foo;
use Cwd;

require Exporter;
@ISA = qw(Exporter);
@EXPORT =
qw(drive hello mergeCoords min max downloadAndDecompress mergeDir downloadRegions);

# This method downloads and decompresses the epigenome files
# It needs the URL for the epigenome, a file that specifically
# Arguments given as ($epigenomeLink, $epigenomeLinkFile, $outputDirectory)
sub drive {
	( my $parentLink, my $inputFile, my $mainDir ) = @_;

	downloadAndDecompress( $parentLink, $inputFile, $mainDir );

	mergeDir($mainDir);
}

# This will eliminate overlaps in a given $inputFile and print out the results into an $outputFile
# If $outputFile is the same as $inputFile, it will be renamed to have an M in it.
# THIS IS THE DEFAULT ACTION AND OTHER SUBROUTINES ARE BASED OFF OF THIS ASSUMPTION!
# Therefore, if $outputFile should be the same as $inputFile if it'll be used in this project.
# The arguments are taken in as ($inputFile, $outputFile)
#
#      It turns this
#
#      ------------     -    ----            -----------      ------     ----------------- -------------
#      ----      ------     -------     - ----       ------   ------------     ------     ------------
#
#      into this
#
#      ---------------- -   -------     - -----------------   ------------------------------------------
sub mergeCoords {

	( my $repeatsInFile, my $repeatsOutFile ) = @_;

	open( IN, "<", $repeatsInFile )
	or die "Cannot open $repeatsInFile: $!\n";

	open( OUT, ">", "$repeatsOutFile" )
	or die "Cannot create $repeatsOutFile: $!\n";

	$_ = <IN>;

	my @line = split;
	while (<IN>) {
		my @temp = split;

		if ( @temp == 0 ) { next; }

		if ( $line[0] eq $temp[0]
			&& isOverlapping( $line[1], $line[2], $temp[1], $temp[2] ) )
		{
			$line[1] = min( $line[1], $temp[1] );
			$line[2] = max( $line[2], $temp[2] );
		}
		else {

			my $result = join( "\t", @line );
			print OUT "$result\n";
			@line = @temp;

		}

	}

	my $result = join( "\t", @line );
	print OUT "$result\n";

	$result = join( "\t", @temp );
	print OUT "$result\n";

	close(IN);
	close(OUT);

}

#
# $inListRef is a list of references to lists
# It looks like this: [[chr, start, end], [chr, start, end], [chr, start, end]]
# This will take the list and merge it.
# Ex.
# chr1	1	10
# chr1	5	15
# chr1	16	20
# chr1	25	30
# chr2	10	20
# chr2	15	25
# chr2	30	40
# =====>
# chr1	1	20
# chr1	25	30
# chr2	10	25
# chr2	30	40
sub mergeCoordsLists {
	( my $inListRef ) = @_;
	my @inList = @$inListRef;

	my @outList = ();
	if ( @inList > 0 ) {
		my @region = @{ $inList[0] };
		my @temp;
		for ( my $i = 1 ; $i <= $#inList ; $i++ ) {
			@temp = @{ $inList[$i] };

			if ( @temp == 0 ) { die "Invalid region\n"; }

			if ( $region[0] eq $temp[0]
				&& isOverlapping( $region[1], $region[2], $temp[1], $temp[2] ) )
			{
				$region[1] = min( $region[1], $temp[1] );
				$region[2] = max( $region[2], $temp[2] );
			}
			else {

				print "@region\n";

				push( @outList, [ $region[0], $region[1], $region[2] ] );
				@region = @temp;

			}
		}
		if ( $region[0] eq $temp[0]
			&& isOverlapping( $region[1], $region[2], $temp[1], $temp[2] ) )
		{
			push( @outList, [ $region[0], $region[1], $region[2] ] );
		}

	}
	return @outList;
}


# returns a list of files in a given directory
# Note: Don't include the '/' at the end
# arguments given as ($parentDirectory)
sub listFiles {
	( my $parentDir ) = @_;
	my $a = `find -f $parentDir`;
	$_ = $a;
	$a =~ s/\/\//\//g;
	return split( "\n", $a );
}

# returns whether or not the epigenomes are overlapping
# based on their starts and ends
# arguments given as ($start1, $end1, $start2, $end2)
sub isOverlapping {
	my ( $s1, $e1, $s2, $e2 ) = @_;

	my $isStartWithin = ( ( $s2 >= $s1 ) && ( $s2 <= $e1 ) );
	my $isEndWithin   = ( ( $e2 >= $s1 ) && ( $e2 <= $e1 ) );
	my $isIncluding   = ( ( $s2 >= $s1 ) && ( $e2 <= $e1 ) );
	my $isIncluded    = ( ( $s1 >= $s2 ) && ( $e1 <= $e2 ) );
	my $isAdjacent = ( $e1 + 1 == $s2 ) || ( $e2 + 1 == $s1 );

	return ( $isStartWithin
		|| $isEndWithin
		|| $isIncluding
		|| $isIncluded
		|| $isAdjacent );
}

# Given a directory, it will perform mergeCoords() on every file
# Arguments given as ($parentDir)
# 
sub mergeDir {
	my ( $fol, $dstDir ) = @_;
	my @files = glob( $fol . '/*.tagAlign' );
	foreach my $file (@files) {
		my @splitFile = split( '/', $file );
		my $nickName = $splitFile[-1];

		if ( $nickName =~ m/.*\d.*tagAlign$/ ) {

			if ( $nickName =~ m/(.*)(\.tagAlign)$/ ) {
				$outFile = "$dstDir/$1M$2";
			}

			print "Merging $file $outFile\n";

			mergeCoords( $file, $outFile );
		}
	}

}

# This runs the unix sort command on epigenomes found
# in a given directory and output the results into another directory.
# Arguments give as (epigenomes directory, output directory)
sub sortEpigenome {
	my ( $epiDir, $sortedEpiDir ) = @_;
	my @markList = glob( $epiDir . "/*.tagAlign" );
	foreach my $mark (@markList) {
		print "\tSorting $mark ...\n";
		my @splitMark = split( '/', $mark );
		my $sortCom =
		"sort -k1,1 -k2,2n $mark > $sortedEpiDir/$splitMark[-1]\n";
		if ( system($sortCom) ) {
			die "Could not execute $sortCom: $!\n";
		}
	}
}

# Simple minimum function
# Arguments give as ($num1, $num2)
sub min {
	my ( $a, $b ) = @_;
	if   ( $a > $b ) { return $b; }
	else             { return $a; }
}

# Simple maximum function
# Arguments give as ($num1, $num2)
sub max {
	my ( $a, $b ) = @_;
	if   ( $a < $b ) { return $b; }
	else             { return $a; }
}

#1;
