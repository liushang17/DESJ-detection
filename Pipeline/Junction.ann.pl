#!/usr/bin/perl -w
use strict;
my $usage = "
Describe:
        Annotation junction with gene names, filtering junctions, creating juntion cell count matrix.
Author:
        liushang\@genomics.cn
Usage:
        perl  $0  <in.gtf>  <Alljunction.filter.list.xls>  <sj.dir>  <Outdir>
\n\n";
die $usage if @ARGV != 4;

my %endid = ();
my %startid = ();
my %endname = ();
my %startname = ();
my %endidtmp = ();
my %startidtmp = ();

##input gtf
my %gtf;
open GTF, $ARGV [0];
while (<GTF>){
       next if /^#/;
       chomp;
       my @t = split /\t/;
       if( $t[2] eq "exon"){
       my ($gid) = $t[-1] =~ /gene_id "([^"]+)"/;
       my ($gname) = $t[-1] =~ /gene_name "([^"]+)"/;
       my $starttmp = $t[4] + 1;
       my $endtmp = $t[3] - 1;
	if((exists $endid{"$t[0]_$endtmp"} and $gid ne $endidtmp{"$t[0]_$endtmp"})){
	$endid{"$t[0]_$endtmp"} .= "-$gid";
	$endname{"$t[0]_$endtmp"} .= "-$gname";
	}
	if(exists $startid{"$t[0]_$starttmp"} and $startidtmp{"$t[0]_$starttmp"} ne $gid){
		$startid{"$t[0]_$starttmp"} .= "-$gid";
		$startname{"$t[0]_$starttmp"} .= "-$gname";
	}
	if(!exists $endid{"$t[0]_$endtmp"} and !exists $startid{"$t[0]_$endtmp"}){
       $endid{"$t[0]_$endtmp"} = $gid;
       $startid{"$t[0]_$starttmp"} = $gid;
       $endname{"$t[0]_$endtmp"} = $gname;
       $startname{"$t[0]_$starttmp"} = $gname;
	}
	$endidtmp{"$t[0]_$endtmp"} = $gid;
	$startidtmp{"$t[0]_$starttmp"} = $gid;
	}
	if($t[2] eq "gene"){
	       my ($gid) = $t[-1] =~ /gene_id "([^"]+)"/;
	       my ($gname) = $t[-1] =~ /gene_name "([^"]+)"/;
		$t[6] = 1 if($t[6] eq '+');
		$t[6] = 2 if($t[6] eq '-');
	       $gtf{$t[0]}{$t[3]}{$t[4]}{$t[6]}{'id'} = $gid;
	       $gtf{$t[0]}{$t[3]}{$t[4]}{$t[6]}{'name'} = $gname;
	}
}
close GTF;

## the annotation of junctions
open AS, $ARGV[1];
open All, ">$ARGV[3]/Alljunction.filter.list.ann.xls" or die "$!";
my @juncinfo;
while (<AS>){
        chomp;
	if ($. == 1){
		print All "junction_id,GeneID,GeneName\n";
		next;
	}
	if($_ =~ /chr/){
        my ($info) =(split /,/)[0];
        my ($chr,$start,$end,$strand) = ("") x 4;
        my @chr = ();
        my @start = ();
        my @end = ();
        if( $info =~ /(chr\w+)_(\d+)_(\d+)_(\d)/g ){
        	$chr = $1;
		$start = $2;
		$end = $3;
		$strand = $4;
	}
        my %over = ();
	my $tmpend = $chr."_".$end;
	my $tmpstart = $chr . "_" . $start;
	if(exists $startid{$tmpstart} and !exists $endid{$tmpend}){
		if(exists $startname{$tmpstart}){
		print All "$_\t$startid{$tmpstart}\t$startname{$tmpstart}\n";
		push @juncinfo, "$_\t$startid{$tmpstart}\t$startname{$tmpstart}\n";
		}else{
			print All "$_\t$startid{$tmpstart}\t-\n";
			push @juncinfo, "$_\t$startid{$tmpstart}\t-\n";
		}
	}elsif(exists $startid{$tmpstart} and exists $endid{$tmpend}){
		if($endid{$tmpend} ne $startid{$tmpstart}){
		if(exists $startname{$tmpstart} and exists $endname{$tmpend}){
		print All "$_\t$endid{$tmpend}:$startid{$tmpstart}\t$endname{$tmpend}:$startname{$tmpstart}\n";
		push @juncinfo, "$_\t$endid{$tmpend}:$startid{$tmpstart}\t$endname{$tmpend}:$startname{$tmpstart}\n";
		}elsif(exists $startname{$tmpstart} and !exists $endname{$tmpend}){
			print All "$_\t$endid{$tmpend}:$startid{$tmpstart}\t-:$startname{$tmpstart}\n";
			push @juncinfo, "$_\t$endid{$tmpend}:$startid{$tmpstart}\t-:$startname{$tmpstart}\n";
		}elsif(!exists $startname{$tmpstart} and exists $endname{$tmpend}){
			print "$_\t$endid{$tmpend}:$startid{$tmpstart}\t$endname{$tmpend}:-\n";
			push @juncinfo,"$_\t$endid{$tmpend}:$startid{$tmpstart}\t$endname{$tmpend}:-\n";
		}else{
			print All "$_\t$endid{$tmpend}:$startid{$tmpstart}\t-+-\n";
			push @juncinfo,"$_\t$endid{$tmpend}:$startid{$tmpstart}\t-+-\n";
		}
		}else{
			if(exists $startname{$tmpstart}){
				print All "$_\t$startid{$tmpstart}\t$startname{$tmpstart}\n";
				push @juncinfo,"$_\t$startid{$tmpstart}\t$startname{$tmpstart}\n";
			}else{
				print All "$_\t$startid{$tmpstart}\t-\n";
				push @juncinfo,"$_\t$startid{$tmpstart}\t-\n";
			}
		}

	}elsif(!exists $startid{$tmpstart} and exists $endid{$tmpend}){
		if(exists $endname{$tmpend}){
		print All "$_\t$endid{$tmpend}\t$endname{$tmpend}\n";
		push @juncinfo,"$_\t$endid{$tmpend}\t$endname{$tmpend}\n";
		}else{
			print All "$_\t$endid{$tmpend}\t-\n";
			push @juncinfo,"$_\t$endid{$tmpend}\t-\n";
		}
	}
	}
}
close AS;

## filter junctions with no or more than 2 genes annotation
open O1,">$ARGV[3]/Alljunction.filter.list.ann.onegene.xls" or die "$!";
my %moregene;
my @junc;
foreach(@juncinfo){
        chomp;
        next if ($_ =~ /junction/);
        my @word = split " ",$_;
        my @sen = split ":",$word[1];
        my $n1 = @sen;
        if($n1 > 1){
                my @sen2 = split "\-",$sen[0];
                my @sen3 = split "\-",$sen[1];
                foreach my $tmp(@sen2){
                        $moregene{$tmp} = 1;
                }
                foreach my $tmp1(@sen3){
                        $moregene{$tmp1} = 1;
                }
        }else{
                my @sen1 = split "\-",$word[1];
                my $n2 = @sen1;
                if($n2 > 1){
                        foreach my $tmp2(@sen1){
                                $moregene{$tmp2} = 1;
                        }
                }else{
                        print  O1 "$_\n";
			my @tmp = split " ",$_;
			push @junc,$tmp[0];
                }
        }
}
close O1;
open O2,">$ARGV[3]/Alljunction.filter.list.ann.moregene.xls" or die "$!";
foreach(keys %moregene){
        print O2 "$_\n";
}
close O2;

print "cell";
foreach(@junc){
        print "\t$_";
}
print "\n";

opendir I1,$ARGV[2] or die "$!";
open O,">$ARGV[3]/merge.count.txt" or die "$!";
foreach(readdir I1){
        chomp;
        next if($_ eq "\." or $_ eq "\.\.");
        my $file = "$ARGV[3]/$_/$_"."SJ.out.tab";
	open I2,$file or die "$!";
        my %sj_num;
        foreach my $line (<I2>){
                chomp($line);
                my @word = split " ",$line;
                my $sj = $word[0]."_".$word[1]."_".$word[2]."_".$word[3];
                my $num = $word[6];
                $sj_num{$sj} = $num;
        }
        close I2;
        print "$_";
        foreach my $junc(@junc){
                if(exists $sj_num{$junc}){
                        print O "\t$sj_num{$junc}";
                }else{
                        print O "\t0";
                }
        }
        print "\n";
}
closedir I1;
close O;
