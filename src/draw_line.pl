#Author : Shujia Huang
#Date   : 2011/02/22
#!/usr/bin/perl
use strict;
use warnings;
use SVG;
use Getopt::Long;
#use Statistics::Basic qw(:all);
use Statistics::Lite qw(:all);
my $MAX_CUTOFF = 80000;
my ( @file1, @file2, $repeat_file, $GC_content, $DMR_file, @common_files );
my ( $window_length, $line_h, $lenged_mark, $title, $help );
my ( $x_manify, $y_manify, $font_manify, $lenged_manify );
my ( $lenged_position );
my ( $depth_cutoff );
GetOptions (
	"a=s"		=> \@file1,
	"b=s"		=> \@file2,
	"dc=i"		=> \$depth_cutoff,
	"GC=s"		=> \$GC_content,
	"D=s"		=> \$DMR_file,
	"r=s"		=> \$repeat_file,
	"f=s"		=> \@common_files,
	"m=s"		=> \$lenged_mark,
	"lp=f"		=> \$lenged_position,
	"lm=f"		=> \$lenged_manify,
	"n=s"		=> \$title,
	"w=i"		=> \$window_length,
	"xm=f"		=> \$x_manify,
	"ym=f"		=> \$y_manify,
	"fm=f"		=> \$font_manify,
	"h=i"		=> \$line_h,
	"help"		=> \$help,
);
sub usage{
	print <<USAGE;
Version : 2.0
Author  : Shujia Huang
Created : 2011/02/22

Last modified : 2012/05/08	Add DMR region. And changed the parameter -GC to -H
Last modified : 2011/02/27
	
		Usage : perl $0 [Options] -a file1 -b file2 -r repeate_file -GC GC_content -f Gene_component_file > output.svg

		Option:

			-a	[str]	bed file. This parameter couldn't be empty.[female or male]
			-b	[str]	bed file. This parameter couldn't be empty.[ male or female]
			-r	[str]	repeate_file. This parameter couldn't be empty.
			-f	[str]	Gene element file. This parameter couldn't be empty.
			-GC	[str]	GC_content file. This parameter couldn't be empty. For heat map.
						
						Input format of -a -b parameters : File:Color1:Color2
						Input format of -r -f parameters : File:Mark:Color
						Input format of -GC   parameters : File:Mark or File:Mark:min:max
			
			-D	[str]	DMR file. For option.	parameter format: File:Mark.
			-m	[str]	lenged mark.
			-lp	[float]	lenged position. [0.4]
			-n	[str]	Title of this figure. [""]
			-xm	[float]	Convergent-divergent x axes. [1]
			-ym	[float]	Convergent-divergent y axes. [1]
			-fm	[float]	Convergent-divergent lenged font size. [1]
			-lm [float] Convergent-divergent lenged size. [1]
			-h	[int]	
			-help	Help

Example:
	perl draw_line.pl -a "A/chr16.coverage.1Mb.xls:rgb(51,153,255):rgb(51,153,255)" -b "J/chr16.coverage.1Mb.xls:rgb(255,0,0):rgb(255,0,0)" -m A:J -r "chromosome_info/new_info/chr16/Repeat.xls.bed:Repeat:rgb(215,215,215)" -f chromosome_info/new_info/chr16/CGI.xls.bed:CGI:"rgb(0,176,80)" -f chromosome_info/new_info/chr16/SNP.xls.bed:SNP:purple -f chromosome_info/new_info/chr16/Gene.xls.bed:Gene:"rgb(255,192,0)" -n "SSC 16" -GC OE.xls.bed:OE:10:60 -ym 0.6 -fm 1.2 -lp 0.75 -lm 1.2 -xm 0.3 -w 1000000 -D chr16.dmr:DMR > chr16.svg
	
USAGE
}
if ( $help or !@file1 or !@file2 or !@common_files or !$GC_content or !$repeat_file){
	usage();
	exit(1);
}
$depth_cutoff  = 40		   if ( !defined $depth_cutoff ); ## Globle value
$x_manify      = 1         if ( !defined $x_manify );
$y_manify      = 1         if ( !defined $y_manify );
$font_manify   = 1 		   if ( !defined $font_manify );
$window_length = 500000    if ( !defined $window_length );
$line_h        = 50        if ( !defined $line_h );
$title         = ""        if ( !defined $title );
$lenged_mark     = "Female:Male" if( !defined $lenged_mark );
$lenged_position = 0.4           if ( !defined $lenged_position);

my @lenged_marks = split ( /:/, $lenged_mark );
my @table1;
my @table2;
my @repeat;
my $colum_num1 = read_file  ( $file1[0], \@table1 );
my $colum_num2 = read_file  ( $file2[0], \@table2 );
my @tissue_inf = shift ( @table1 ); 
shift ( @table2 );
die "Tissue unmatch!!!!\n" if ( $colum_num1 != $colum_num2 );

my $max_col    = ( $colum_num1 >= $colum_num2 ) ? $colum_num1 : $colum_num2;
my $max_row    = ( $#table1 >= $#table2 ) ? $#table1: $#table2;
my $x_length   = $x_manify * $window_length * $max_row / 100000;
my $y_length = 480 * $y_manify;
my $x_flank  = 50;
my $y_flank  = 50;
my $x_start  = 80;
my $y_start  = 60 + 1.2 * ($line_h + 10);
my $chr_width  = 20 * $y_manify;
my $d_band     = 10;
my $band_width = 10;
my $h          = 50  * $y_manify;
my $fig_x_size = $x_start + $x_length  + $x_flank; 
my $fig_y_size = $y_start + $chr_width + $h + $d_band + $band_width + $y_flank;
my $distance   = ( $lenged_position <= 1) ? $lenged_position * $x_length : ( $x_length + 10 );
my $svg =  SVG->new('width', $fig_x_size,'height', $fig_y_size );
my ( $max_cg, $min_cg ) = draw_heatmap ( $x_start, $y_start, $chr_width, $x_length, ($x_length / $max_row), $GC_content, 0, \$svg );

my ( $sub_file, $mark, $color ) = split ( /:/, $repeat_file );
read_file ( $sub_file, \@repeat );
my $max = get_extremum( \@repeat, 3 );
my $element_num = @common_files;
my $d           = $h / ($element_num + 1);
my $font_size   = int( 0.8 * $d - 1 ) + 1;
my $y_text      = $y_start + $chr_width + $d;
#my $text_color  = ( $color eq "none" ) ? "black" : $color;
my $text_color  = "black";
$svg->text( 'x' => $x_start - 0.6 * $font_size * length($mark) - 1, 'y' => $y_text, '-cdata', $mark, 'fill', $text_color, 'font-family','Arial','font-size', $font_size );
draw_histogram( $x_start, $y_start + $chr_width + $h, $x_length / $max_row, 0.8 * $h, $max, $window_length, $color, \@repeat, \$svg );

for ( my $i = 0; $i < $element_num; ++$i ){
	( $sub_file, $mark, $color ) = split ( /:/, $common_files[$i] );
	my @table;
	my $y   = $y_start + $chr_width + $h;
	$y_text = $y_start + $chr_width + ($i + 1 + 1 ) * $d;

	read_file ( $sub_file, \@table );
	my $max = get_extremum( \@table );
	$svg->text( 'x' => $x_start - 0.6 * $font_size * length($mark) - 5, 'y' => $y_text, '-cdata', $mark, 'fill', $color, 'font-family','Arial','font-size', $font_size );
	draw_line( $x_start, $y, $x_length / $max_row, $max, $window_length, 0.8 * $h, \@table, -1, $color, $color, 0, 0, \$svg);
}
my $max2  = get_file_maxs( @file1, @file2 );
$max2     = $depth_cutoff;
#my $flag = draw_up( $x_start, $y_start, $x_length, $max_row, $max2, $window_length, $line_h, $chr_width + $h, $distance, $lenged_marks[0], 1, \@file1, \$svg, 1 );
#draw_up( $x_start, $y_start, $x_length, $max_row, $max2, $window_length, $line_h, $chr_width + $h, $distance, $lenged_marks[1], 0, \@file2, \$svg, $flag );
my $color1 = draw_up( $x_start, $y_start, $x_length, $max_row, $max2, $window_length, $line_h, 1, \@file1, \$svg );
my $color2 = draw_up( $x_start, $y_start, $x_length, $max_row, $max2, $window_length, $line_h, 0, \@file2, \$svg );

## Draw DMR line
draw_heatmap ( $x_start, $y_start + $chr_width + $h + $d_band * 0.5, 0.75 * $chr_width, $x_length, ($x_length / $max_row), $DMR_file, 1, \$svg ) if ( defined $DMR_file );
#draw_heatmap ( $x_start, $y_start + $chr_width + $h + $d_band * 0.5, 0.75 * $chr_width, $x_length, ($x_length / $max_row), $DMR_file, 0, \$svg ) if ( defined $DMR_file );
draw_color_band( $x_start, $y_start + $chr_width + $h + $d_band + $chr_width, $band_width, $max_cg, $min_cg, \$svg );
draw_ruler ( $x_start,     $y_start + $chr_width + $h + $d_band + $chr_width + $band_width + 25, 5 * ($x_length / $max_row), "black", "5M", \$svg );

$font_size = ( $font_manify < 1 ) ? 10 * $font_manify : 10;
draw_lenged( $x_start + $distance, $y_start + $chr_width + $h + $d_band + $chr_width, $color1, $font_size, $font_manify, $lenged_marks[0], \$svg );
draw_lenged( $x_start + $distance, $y_start + $chr_width + $h + $d_band + $chr_width + 1.5 * $font_size, $color2, $font_size, $font_manify, $lenged_marks[1], \$svg );

draw_title ( $x_start, $y_start -($line_h + 10) - 20, $x_length, $title, \$svg );
print $svg->xmlify();
#################################################################################################################################
sub get_file_maxs {
	my ( @files ) = @_;
	my @max;
	for ( my $i = 0; $i < @files; ++$i ){
		my @table;	
		my ( $sub_file ) = (split ( /:/, $files[$i] ))[0];
		read_file ( $sub_file, \@table );
		shift ( @table );
		my @maxs = get_max( \@table );
    	push ( @max, max( @maxs ) );
	}
	my $max = max( @max );
	
	return $max;
}
sub draw_up {
#	my ( $x_start, $y_start, $x_length, $max_row, $max, $window_length, $line_h, $h, $distance, $mark, $is_mark, $file, $svg, $flag ) = @_;
	my ( $x_start, $y_start, $x_length, $max_row, $max, $window_length, $line_h, $is_mark, $file, $svg ) = @_;
	my $file_num  = @$file;
	#my $font_size = ( $font_manify < 1 ) ? 10 * $font_manify : 10;
	my ( $x, $y, @table );
	my ( $sub_file, $color1, $color2, $breed );
	for ( my $i = 0; $i < $file_num; ++$i ){
		( $sub_file, $color1, $color2, $breed ) = split ( /:/, $$file[$i] );
		$breed  = ( !defined $breed  ) ? "" : $breed;
		$color1 = ( !defined $color1 ) ? "red" : $color1;
		$color2 = ( !defined $color2 ) ? "blue" : $color2;
		@table  = ();
		read_file ( $sub_file, \@table );
		shift ( @table );
		draw_line ( $x_start, $y_start - 10 , $x_length / $max_row, $max, $window_length, $line_h, \@table, -1, $color1, $color2, 1, $is_mark, $svg);
		$is_mark = 0;
	}
=head
	if ( $flag ){
		draw_lenged( $x_start + $distance, $y_start + 10 + $h + 1.5 * $font_size, $color1, $font_size, $font_manify, "Fat", $svg );
		draw_lenged( $x_start + $distance, $y_start + 10 + $h, $color2, $font_size, $font_manify, "Muscle", $svg );
		#draw_lenged( $x_start + $distance, $y_start - ( $line_h + 10 ), $color1, $font_size, $font_manify, "Fat", $svg );
		#draw_lenged( $x_start + $distance, $y_start - ( $line_h + 10 ) + 1.5 * $font_size, $color2, $font_size, $font_manify, "Muscle", $svg );
	}
	return 0;
=cut
	return ( $color1 );
}

sub read_file {
	my ( $file, $table ) = @_;
	($file) = split ( /:/, $file );
	my @tmp;
	my $colum_num;
	open ( IN, $file ) or die "Cannot open : $file\n";
	while ( <IN> ) {
		chomp;
		@tmp = split ( /\s+/, $_ );
		next if ( $tmp[1] !~ m/^\d+$/ ); 
		push ( @$table, [ @tmp ] );
		$colum_num = @tmp;
	}
	close ( IN );
	return ($colum_num - 3);
}
sub draw_heatmap {
	my ( $x_start, $y_start, $height, $width, $delta, $file_inf, $no_grade, $svg ) = @_;
	my @Inf = split( /:/, $file_inf );
	die "[Format ERROR]The Format of -GC or -D parameter ERROR. Not like this: $file_inf\n" if( @Inf < 2 );
	my ($file, $mark, $min, $max ) = @Inf;
	my @table;
	read_file( $file, \@table );

	if ( $no_grade ) { 
		$max = 1; 
		$min = 1; 
	} else {
		$max = get_extremum( \@table, 3 ) if ( !defined $max );
	    $min = get_min( \@table )         if ( !defined $min );
	}

	my $grade       = ( $min < $max ) ? (int( $max - $min ) + 1) : 1;
	my $font_size   = int( 0.6 * $height ) + 1;
	my $font_col    = ( $no_grade ) ? "black" : "gray";
	$$svg->text( 'x', $x_start - 0.8 * $font_size * length($mark) -5 , 'y', $y_start + 0.8 * $height, '-cdata', $mark, 'fill', $font_col, 'font-family','Arial','font-size', $font_size );
#	$$svg->rect( 'x', $x_start, 'y', $y_start, 'width', $width, 'height', $height, 'stroke', "rgb(215,215,215)", 'fill', "none", 'stroke-width', 2 );
	for ( my $i = 0; $i < @table; ++$i ){
		my $x   = $x_start + $delta * ( int( $table[$i][1] / $window_length ) );
		my $y   = $y_start;
		my $color;
		my $col;
		if ( $no_grade ) {
			$col = 0; # Black
		} else {
			$col = int ( 255 * ( 1 - (int($table[$i][3] - $min )) / $grade ));
		}
		$color = "rgb( $col, $col, $col )";
#my $w  = ( $no_grade ) ? ( $table[$i][2] - $table[$i][1] + 1 ) * $delta / $window_length : $delta; 
my $w  = ( $no_grade ) ? 0.1: $delta; 
		$$svg->rect( 'x', $x, 'y', $y, 'width', $w, 'height', $height, 'fill', $color );

#		if ( $i == 7 ){
#			$$svg->rect( 'x', $x - $delta * 7, 'y', $y, 'width', $delta * 7, 'height', $height, 'fill', "none", 'stroke', "red" );
#		}
#		if ( $i == @table - 7 ){
#			$$svg->rect( 'x', $x, 'y', $y, 'width', $delta * 7, 'height', $height, 'fill', "none", 'stroke', "red" );
#		}
	}
	return ( $max, $min );
}
sub draw_color_band {
	my ( $x_start, $y_start, $band_width, $max, $min, $svg ) = @_;
	my $length = 60;
	my $grade  = ( $min < $max ) ? (int( $max - $min ) + 1) : 1; $grade = int ( $max ) + 1;
	my $delta  = $length / $grade;
	my $font_size = ( $font_manify < 1) ? 10 : 10 * $font_manify;
	$max       = sprintf "%.1f", $max; $max .= "%";
	$min	   = sprintf "%.1f", $min; $min .= "%";
	###
	#$min = "30%";
	#$max = "60%";
	$$svg->rect( 'x' => $x_start,'y' => $y_start, 'width', $length, 'height', $band_width, 'stroke', 'black', 'fill', "none" );
	$$svg->line ( 'x1'=> $x_start, 'y1'=> $y_start + $band_width, 'x2' => $x_start, 'y2' => $y_start + $band_width + 3, 'stroke', 'gray');
	$$svg->line ( 'x1'=> $x_start + $length, 'y1'=> $y_start + $band_width, 'x2' => $x_start + $length, 'y2' => $y_start + $band_width + 3, 'stroke', 'gray');
	$$svg->text( 'x', $x_start - 0.25 * length( $min ) * $font_size, 'y', $y_start + $band_width + $font_size + 3, '-cdata', $min, 'fill', "black",'font-family','Arial','font-size', $font_size );
	$$svg->text( 'x', $x_start + $length - 0.25 * length( $max ) * $font_size, 'y', $y_start + $band_width + $font_size + 3, '-cdata', $max, 'fill', "black",'font-family','Arial','font-size', $font_size );
	###
	for ( my $i = 0; $i < $grade; ++$i ){
		my $x   = $x_start + $delta * $i;
		my $y   = $y_start;
		my $col = int ( 255 * ( 1 - ( $i + 1 ) / $grade ) );
		my $color = "rgb( $col, $col, $col )";
		$$svg->rect( 'x' => $x, 'y' => $y, 'width', $delta, 'height', $band_width, 'fill', $color, 'stroke', $color );
	}
    ###
} #
sub draw_histogram {
	my ( $x_start, $y_start, $delta, $h, $max, $window_length, $color, $table, $svg ) = @_;
	for ( my $i = 0; $i <@$table; ++$i ){
		my $x     = $x_start + $delta * ( int( $$table[$i][1] / $window_length ) );
		my $new_h = $h * $$table[$i][3] / $max;
		$$svg->rect( 'x' => $x, 'y' => $y_start - $new_h, 'width', $delta, 'height', $new_h, 'fill', $color, 'stroke', 'black','stroke-width', 0.3 );
	}
}
sub get_max {
	my ( $table ) = @_;
	my @max;

	my $column = @{$$table[0]};
	for ( my $i = 3; $i < $column; ++$i ){
		$max[$i-3] = get_extremum( $table, $i );
	}
	return ( @max );
}
sub draw_line {
	my ( $x_start, $y_start, $delta, $max, $window_length, $sub_h, $table, $upDown, $color1, $color2, $is_small, $is_mark, $svg ) = @_; ##upDown 0 or 1
	my @init_coor;
	my $flag = 1;
	my @y1   = @{$$table[0]};
	my $pos1 = int( $$table[0][1] / $window_length );
	my $x1   = $x_start + $delta * ( $pos1 + 0.5 );

#	$max = ( $max > $depth_cutoff ) ? $depth_cutoff : $depth_cutoff if ( $is_small );
	for ( my $i = 0; $i < @$table; ++$i ) {
		my $pos2 = int($$table[$i][1] / $window_length);
		my $x2   = $x_start + $delta * ( $pos2 + 0.5 );
		$x1      = ( $pos2 > $pos1 + 1 ) ? $x2 : $x1;
		@y1      = ( $pos2 > $pos1 + 1 ) ? @{$$table[$i]} : @y1;

		my $color;
		for ( my $k = @{$$table[$i]} - 1, my $j = 3; $j < @{$$table[$i]}; ++$j, --$k ) {
			my $over        = ( $$table[$i][$j] >= $MAX_CUTOFF ) ? 1 : 0;
			my $new_over    = ( $$table[$i][$j] > $depth_cutoff ) ? 0 : 0;
			#$y1[$j]         = ( $y1[$j] >= $MAX_CUTOFF ) ? $max : $y1[$j];			
			#$$table[$i][$j] = ( $$table[$i][$j] >= $MAX_CUTOFF ) ? $max : $$table[$i][$j];			
			$y1[$j]         = ( $y1[$j] >= $depth_cutoff ) ? $depth_cutoff : $y1[$j] if ( $is_small );
			$$table[$i][$j] = ( $$table[$i][$j] >= $depth_cutoff ) ? $depth_cutoff : $$table[$i][$j] if ( $is_small );
			my $tmp_h1= $y1[$j] * $sub_h / $max;
			my $tmp_h2= $$table[$i][$j] * $sub_h / $max;
			$tmp_h1  = ( $tmp_h1 > $sub_h ) ? $sub_h : $tmp_h1;
			$tmp_h2  = ( $tmp_h2 > $sub_h ) ? $sub_h : $tmp_h2;
			my $y1   = $y_start + $upDown * $tmp_h1;
			my $y2   = $y_start + $upDown * $tmp_h2;
			my $size = ( $is_small ) ? 0.2 : 1;
			$size    = ( $new_over ) ? 0.8 : $size; 
			$color = ( $k > 4 ) ? $color1 : $color2;

			$$svg->circle( 'cx', $x2, 'cy', $y2, 'r', $size, 'fill', $color );
			if ( !$new_over || !$is_small ){ ## $is_small means draw the 60 tissues line!!
				$$svg->line( 'x1' => $x1, 'y1' => $y1, 'x2' => $x2, 'y2' => $y2, 'stroke', $color, 'stroke-width', $size );
			} 			
=head	
			if ( $over ){
				### If bigger than MAX_CUTOFF, draw 3 little circle above this point!
				for ( my $m = 0; $m < 3; ++$m ){
					my $y = $y2 + $upDown * ( $m + 1 ) * $sub_h / 10;
					$$svg->circle( 'cx', $x2, 'cy', $y, 'r', 0.5, 'fill', "gray" );
				}
			}
=cut

			if ( $flag ) {
				$y2 = $y_start + $upDown * ( ( $j - 3 ) * ( $sub_h + 10 ) + $sub_h );
				$$svg->line( 'x1', $x_start, 'y1', $y2 - $upDown * $sub_h, 'x2', $x_start, 'y2', $y2 + $upDown * $delta,'stroke', "rgb(215,215,215)" );
				###
				$$svg->line( 'x1', $x_start, 'y1', $y2 - $upDown * $sub_h, 'x2', $x_start - 2, 'y2', $y2 - $upDown * $sub_h, 'stroke', "rgb(215,215,215)" );
				$$svg->line( 'x1', $x_start, 'y1', $y2 + $upDown * $delta, 'x2', $x_start - 2, 'y2', $y2 + $upDown * $delta, 'stroke', "rgb(215,215,215)" );
				if ( $is_mark ) {
					my $font_size = ( $font_manify < 0 ) ? 10 : 10 * $font_manify;
					my $mark = ( $max > 1000 ) ? (sprintf "%.1e", $max ) : $max;
					$$svg->text( 'x', $x_start - 0.5 * length( $mark ) * $font_size - 5, 'y', $y2 + $upDown * $delta + 0.4 * $font_size, '-cdata', $mark, 'fill', "black",'font-family','Arial','font-size', $font_size );
					$$svg->text( 'x', $x_start - 0.5 * length( '0' ) * $font_size - 5, 'y', $y2 - $upDown * $sub_h + 0.4 * $font_size, '-cdata', '0', 'fill', "black",'font-family','Arial','font-size', $font_size );
					$is_mark = 0;
				}
				###
				push ( @init_coor, [ $x_start, $y2, $k ]);
				$flag = 0;
			}
#			$flag = 0;
		}
		@y1   = @{$$table[$i]};	
		$x1   = $x2;
		$pos1 = $pos2;
	}
	return ( @init_coor );
}
sub draw_tissue_name {
	my ( $tissue_inf, $tissue_coor, $line_h, $upDown, $svg ) = @_;

	my $font_size = int (0.5 * $line_h) + 1;
	for ( my $i = 0; $i < @$tissue_coor; ++$i ) {
		my $x = $$tissue_coor[$i][0] - 2.8 * $font_size;
		my $y = $$tissue_coor[$i][1] -$upDown * 0.5 * $line_h;
		my $name = $$tissue_inf[0][ $$tissue_coor[$i][2] ];
		$$svg->text( 'x' => $x, 'y' => $y, '-cdata', $name, 'fill', 'black', 'font-family','Arial','font-size', $font_size);
	}
}
sub draw_lenged_old {
	my ( $x_start, $y_start, $color1, $color2, $manify, $svg ) = @_;
	my $line_length = ( $manify > 1 ) ? 30 : 30 * $manify;
	my $crevice     = ( $manify > 1 ) ? 10 : 10 * $manify;
		
	$$svg->line( 'x1' => $x_start, 'y1' => $y_start, 'x2' => $x_start + $line_length, 'y2' => $y_start, 'stroke', $color1 , 'stroke-width', 2 );
	$$svg->line( 'x1' => $x_start + $line_length + $crevice, 'y1' => $y_start, 'x2' => $x_start + 2 * $line_length + $crevice, 'y2' => $y_start, 'stroke', $color2, 'stroke-width', 2 );

	return ( $x_start + 2 * ( $line_length + $crevice ), $y_start );
}

sub draw_lenged {
    my ( $x_start, $y_start, $color, $font_size, $manify, $text, $svg ) = @_;
	my $line_length = ( $manify > 1 ) ? 30 : 30 * $manify;
	my $crevice     = ( $manify > 1 ) ? 10 : 10 * $manify;

	$$svg->line('x1' => $x_start, 'y1' => $y_start, 'x2' => $x_start + $line_length, 'y2' => $y_start, 'stroke', $color, 'stroke-width', 2);
	$$svg->text( 'x', $x_start + $line_length + $crevice, 'y', $y_start + 0.4 * $font_size, '-cdata', $text, 'fill', 'black', 'font-family','Arial','font-size',$font_size );
}

sub draw_ruler {
    my ( $x_start, $y_start, $line_length, $color, $text, $svg ) = @_;
	my $font_size = ( $font_manify < 1 ) ? 10 * $font_manify : 10;
    my $crevice   = ( $font_manify > 1 ) ? 10 : 10 * $font_manify;

    $$svg->line('x1' => $x_start, 'y1' => $y_start, 'x2' => $x_start + $line_length, 'y2' => $y_start, 'stroke', $color, 'stroke-width', 2);
	$$svg->line('x1' => $x_start, 'y1' => $y_start, 'x2' => $x_start, 'y2' => $y_start - 3, 'stroke', $color, 'stroke-width', 2);
	$$svg->line('x1' => $x_start + $line_length, 'y1' => $y_start, 'x2' => $x_start + $line_length, 'y2' => $y_start - 3, 'stroke', $color, 'stroke-width', 2);

    $$svg->text( 'x', $x_start + $line_length + $crevice, 'y', $y_start + 0.4 * $font_size, '-cdata', $text, 'fill', 'black', 'font-family','Arial','font-size',$font_size );
}

sub draw_title {
	my ( $x_start, $y_start, $x_length, $text, $svg ) = @_;
	my $text_length = length( $text );
#	my $font_size   = 16;
	my $font_size = ( $font_manify < 1 ) ? 10 * $font_manify : 10;
	my $x           = $x_start + 0.5 * ($x_length - 0.5 * $text_length * $font_size);
	my $y           = $y_start - 0.5 * $font_size;
	
	$$svg->text( 'x' => $x, 'y' => $y, '-cdata', $text, 'fill', 'black', 'font-family','Arial','font-size', $font_size );
}
sub get_extremum {
	my ( $table, $index ) = @_;
	$index = 3 if ( !defined $index);
	die "Table element ERROR\n" if( @{$$table[0]} < 4 );
	my $max = $$table[0][$index];
	for ( my $i = 1; $i < @$table; ++$i ){
		### Get max value and ensure max value is lower than MAX_CUTOFF
		$max = ( $max < $$table[$i][$index] && $$table[$i][$index] < $MAX_CUTOFF ) ? $$table[$i][$index] : $max;
	}
	return $max;
}
sub get_min {
    my ( $table ) = @_;
    die "Table element ERROR\n" if( @{$$table[0]} < 4 );
    my $min = $$table[0][3];
    for ( my $i = 1; $i < @$table; ++$i ){
        $min = ( $min > $$table[$i][3] ) ? $$table[$i][3] : $min;
    }
    return $min;
}
