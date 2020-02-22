#!/usr/bin/perl -w
use strict;
use warnings;

use FindBin qw($Bin);
use File::Basename;
use File::Spec::Functions qw(catdir);
use Cwd;
use lib catdir(dirname($Bin), '/lib');
use Getopt::Long;
use Capture::Tiny ':all';
use IO::Compress::Gzip qw(gzip $GzipError);
use Log::Log4perl;
use Cwd;


die ("Usage: FR-Hit_parse.pl <Output_FR-Hit> <ContigFile> <prefix> <OUTFILE>\n") unless (scalar(@ARGV) == 4);

use constant WAIT_TIME => '60' ; # 60 seconds = 1 minute
use constant NO_CHANGE_IN_FILE_FOR_X_SECONDS => '30';

# perl $Bin/FR-Hit_link_light_v2.pl $frhitFile $CDHITFILE $sampleId $OUTFILE="ForPhrap2.fasta"
# noVector ==> phrap.cdhit1
# file == frhit
my $file = shift @ARGV;
# fna == cdhit
my $fna = shift @ARGV;
my $prefix = shift @ARGV;
my $output_file = shift @ARGV;
my $dirname = dirname($file);
my $pwd;
my @seqs=();
my %links=();
my $linkCount = 0;
my $linkThreshold = 5;
my %linkNum=();

open (IN, "<$file") or die ("Couldn't open file: $file\n");
open (FNA, "<$fna") or die ("Couldn't open file: $fna\n");

my %len=();
my %refs=();
my $name="";
# 读取去重后的phrap组装序列
while (my $k = <FNA>){
    chomp($k);
    if ($k =~ /^>(\S+)\s*$/){
        print "Name is >$name< and size is $len{$name}\n";
        $name=$1;
    }else{
        $refs{$name}.=$k;
        $len{$name}=length($refs{$name});
    }
}


my $old="";
my $n="";
# 读取noVector序列映射到phrap序列的映射关系
# 得到FW和RC映射到不同phrap序列的连接方式
while (my $l = <IN>){
    chomp ($l);
    my @array = split /\t/, $l;
    $n = $array[0];
    $n =~ s/\_\d$//;

    $array[7] =~ s/\%$//;
    $array[1] =~ s/nt$//;
    # 相似度小于95%跳过该序列
    next if ($array[7] < 95);

    if ($n ne $old){
        if ($old ne ""){
            &process_array(\@seqs, $n);
        }
        $old = $n;
        @seqs=();
    }
    push (@seqs, [@array]);
}
close IN;
&process_array(\@seqs, $n);

my %contigs=();
my %used=();
my %skip=();
my $counter=1;
my $count2=1;

#die ("End of Test\n");

foreach my $l (keys(%links)){
    foreach my $t (keys(%{$links{$l}})){
        if ($links{$l}{$t} >= $linkThreshold){
            $linkCount++;
            $linkNum{$l}{$t} = $linkCount;
        }
    }
}

if ($linkCount > 39){
    $linkThreshold = 10;
    $linkCount = 0;
    foreach my $l (keys(%links)){
        foreach my $t (keys(%{$links{$l}})){
            if ($links{$l}{$t} >= $linkThreshold){
                $linkCount++;
                $linkNum{$l}{$t} = $linkCount;
            }
        }
    }
    if ($linkCount > 39){
        $linkThreshold = 20;
        $linkCount = 0;
        foreach my $l (keys(%links)){
            foreach my $t (keys(%{$links{$l}})){
                if ($links{$l}{$t} >= $linkThreshold){
                    $linkCount++;
                    $linkNum{$l}{$t} = $linkCount;
                }
            }
        }
        if ($linkCount > 39){
            $linkThreshold = 40;
            $linkCount = 0;
            foreach my $l (keys(%links)){
                foreach my $t (keys(%{$links{$l}})){
                    if ($links{$l}{$t} >= $linkThreshold){
                        $linkCount++;
                        $linkNum{$l}{$t} = $linkCount;
                    }
                }
            }
            if ($linkCount > 39){
                $linkThreshold = 60;
                $linkCount = 0;
                foreach my $l (keys(%links)){
                    foreach my $t (keys(%{$links{$l}})){
                        if ($links{$l}{$t} >= $linkThreshold){
                            $linkCount++;
                            $linkNum{$l}{$t} = $linkCount;
                        }
                    }
                }
                if ($linkCount > 39){
                    $linkThreshold = 80;
                    $linkCount = 0;
                    foreach my $l (keys(%links)){
                        foreach my $t (keys(%{$links{$l}})){
                            if ($links{$l}{$t} >= $linkThreshold){
                                $linkCount++;
                                $linkNum{$l}{$t} = $linkCount;
                            }
                        }
                    }
                    if ($linkCount > 39){
                        $linkThreshold = 100;
                        $linkCount = 0;
                        foreach my $l (keys(%links)){
                            foreach my $t (keys(%{$links{$l}})){
                                if ($links{$l}{$t} >= $linkThreshold){
                                    $linkCount++;
                                    $linkNum{$l}{$t} = $linkCount;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

$pwd = cwd();
print "Current_Dir: $pwd\n";
my $linkDir = "$pwd/linkFiles";
if (not(-d $linkDir)) {
    system("mkdir $linkDir");
}
my $toMap_file = "$linkDir/$prefix-Link_ToMap.fna";
my $linkMap_file = "$linkDir/$prefix-Link_Map.txt";
my $toGet_Seq_file = "$linkDir/$prefix-Link_ToGet_Seqs.txt";
my $linkOut_file = "$linkDir/$prefix-Link_OutCdHit.fna";
my $noChimera_file = "$linkDir/$prefix-Link_NoChimera.fna";

foreach my $l (keys(%links)){
    foreach my $t (keys(%{$links{$l}})){
        if ($links{$l}{$t} >= $linkThreshold){
            print "Evaluating link between $l and $t that has $links{$l}{$t}\n";
            print "Link number $linkNum{$l}{$t} has a count of $links{$l}{$t} reads. A total of $linkCount links are being evaluated and the threshold for evaluation is $linkThreshold reads\n";
            my @gens = split /\!/, $l;
            my @pos = split /\!/, $t;
            my $sub1="";
            my $sub2="";
            my $return_file="";

            $sub1 = substr $refs{$gens[0]}, 0, 100 if ($pos[0] =~/start/);
            $sub1 = substr $refs{$gens[0]}, -100 if ($pos[0] =~/stop/);
            $sub2 = substr $refs{$gens[1]}, 0, 100 if ($pos[1] =~/start/);
            $sub2 = substr $refs{$gens[1]}, -100 if ($pos[1] =~/stop/);
            $used{"$gens[0].$counter"}=$refs{$gens[0]};
            $used{"$gens[1].$counter"}=$refs{$gens[1]};
            $skip{$gens[0]}=1;
            $skip{$gens[1]}=1;

            open (OUT, ">$toMap_file") or die ("Couldn't open file $toMap_file\n");
            print OUT ">$pos[0]\_$gens[0]\n$sub1\n>$pos[1]\_$gens[1]\n$sub2\n";
            close OUT;
            # 以上100个碱基Link序列映射到去接头引物的Fasta文件
            system("fr-hit -d $toMap_file -a $dirname/$prefix.fasta.clean -o $linkMap_file -g 1 -q 50 -c 95");
            # 获得映射错误的序列名
            system("perl /home/ubuntu/Parfums/sub/FR-Hit_cleanLink.pl $linkMap_file $toGet_Seq_file");
            # 根据映射错误的序列名得到对应的序列，追加到Link序列文件中
            system("perl /home/ubuntu/Parfums/sub/GetSequences.pl $toGet_Seq_file $dirname/$prefix.fasta.clean >> $toMap_file");
            # 去重
            system("cd-hit-est -i $toMap_file -o $linkOut_file -G 0 -aS 0.99 -g 1 -r 1 -c 0.9");
            my $len=60;
            # 组装去重后的序列
            system("phrap -minmatch 10 -maxmatch 30 -bandwidth 0 -minscore 15 $linkOut_file");
            # 将组装去重后的序列映射到去接头引物的Fasta文件
            system("fr-hit -d $linkOut_file.contigs -o $linkMap_file -a $dirname/$prefix.fasta.clean -m 30");
            # 去嵌合体
            system("perl /home/ubuntu/Parfums/sub/FR-Hit_cleanChimera.pl $linkMap_file $linkOut_file.contigs > $noChimera_file");
            my @contig =`cat $noChimera_file`;
            my $na = "";
            my $contigsRound=0;
            my $notPass="";
            # 如果去嵌合体后的序列长度小于$len,丢弃该序列
            foreach my $s (@contig){
                chomp($s);
                if ($s =~ /^>/){
                    if ($na ne ""){
                        if (length($contigs{$na}) < $len){
                            $notPass.=">$na\n$contigs{$na}\n";
                            delete $contigs{$na};
                        }else{
                            $contigsRound++;
                        }
                        #print "Contig is \n$na\n$contigs{$na}\n";
                    }
                    $na=$s;
                    $na.=".$counter";
                    $counter++;
                }else{
                    $contigs{$na}.=$s;
                }
            }
            if (length($contigs{$na}) < $len){
                $notPass.=">$na\n$contigs{$na}\n";
                delete $contigs{$na};
            }else{
                $contigsRound++;
            }
            #print "Contig is \n$na\n$contigs{$na}\n";

            print "Testing link No $count2 of $linkCount total links\n";
            $count2++;


            #Let's do some cleaning

            #die "No Contigs pass the threshold, len was $len, contigs are:\n$notPass\n" unless ($contigsRound > 0);
            #system("rm -r Link* SGE");
        }
    }
}

open (OUT, ">$output_file") or die ("Couldn't open file $output_file\n");

foreach my $o (keys(%refs)){
    # 输出没有被修改的原序列
    print OUT ">$o\n$refs{$o}\n" unless ($skip{$o});
}
foreach my $p (keys(%used)){
    # 输出被修改过且序列长度大于30的序列
    print OUT ">$p\n$used{$p}\n" unless length($used{$p}) < 30;
}
foreach my $u (keys(%contigs)){
    # 输出组装后的序列
    print OUT "$u\n$contigs{$u}\n" if ($contigs{$u} && length($contigs{$u}) > 30);
}
close OUT;


########### subroutines #######

sub process_array{
    my @arr = @{$_[0]};
    my $name = $_[1];
    my %seqs=();
    my %gens=();
    # $seqs{$arr[$i][8]}{$arr[$i][0]} == {被映射序列, 映射序列}
    for (my $i=0; $i<@arr; $i++){
        # 以字典为数组名，添加映射关系
        push (@{$seqs{$arr[$i][8]}{$arr[$i][0]}}, [@{$arr[$i]}]);
        # {映射序列，被映射序列} 映射次数计数+1
        $gens{$arr[$i][0]}{$arr[$i][8]}++;
    }
    # 如果映射序列 的存在个数不等于2，跳出
    # 即FW和RC端都有映射关系
    my @ke = keys(%gens);
    return unless scalar(@ke ==2);
    # geA[0], geB[0]两条不同的被映射序列
    # 如果被映射序列为同一条，跳出
    my @geA = keys (%{$gens{$ke[0]}});
    my @geB = keys (%{$gens{$ke[1]}});
    return unless (scalar(@geA)==1 && scalar(@geB)==1 && $geA[0] ne $geB[0]);


    my $switch1=0;
    my $l=0;

    # gens字典必须是有序排列的，取出时对应前后端?
    # Python使用OrderedDict
    my @a1=@{${$seqs{$geA[0]}{$ke[0]}}[0]};
    my @a2=@{${$seqs{$geB[0]}{$ke[1]}}[0]};

    die ("No existe len de >$geA[0]<\n") unless ($len{$geA[0]});
    die ("No existe len de >$geB[0]<\n") unless ($len{$geB[0]});


    my $link="";
    if ($a1[6] eq "+"){
        # 正向映射，结束位置大于被映射序列的长度-100
        if($a1[-1] > ($len{$geA[0]}-100)){
            $link="stop_1";
        }else{
            $link="weird1";
        }
    }else{
        # 逆向映射，映射结束位置小于100
        if ($a1[-2] < 100){
            $link="start_1";
        }else{
            $link="weird_1";
        }
    }
    if ($a2[6] eq "+"){
        if($a2[-1] > ($len{$geB[0]}-100)){
            $link.="!stop_2";
        }else{
            $link.="!weird2";
        }
    }else{
        if ($a2[-2] < 100){
            $link.="!start_2";
        }else{
            $link.="!weird_2";
        }
    }

    #print "@a1\n@a2\n$link\nWhere\n" if ($geA[0] eq "Clean_seq_964_8" && $geB[1] eq "Clean_seq_169_11" && $link =~ /stop_1-start_2/);
    # 仅处理两条序列都映射到被映射序列两边的情况
    if ($link !~ /weird/){
        # 情况一：映射都为同向
        if ($link eq "start_1!start_2" or $link eq "stop_1!stop_2" ){
            if ($links{"$geB[0]!$geA[0]"}{$link}){
                $links{"$geB[0]!$geA[0]"}{$link}++;
            }else{
                $links{"$geA[0]!$geB[0]"}{$link}++;
            }
        # 情况二：FW为逆向，RC为正向
        }elsif ($link eq "start_1!stop_2"){
            if ($links{"$geB[0]!$geA[0]"}{"stop_1!start_2"}){
                $links{"$geB[0]!$geA[0]"}{"stop_1!start_2"}++;
            }else{
                $links{"$geA[0]!$geB[0]"}{$link}++;
            }
        # 情况三：FW为正向，RC为逆向
        }elsif ($link eq "stop_1!start_2"){
            if ($links{"$geB[0]!$geA[0]"}{"start_1!stop_2"}){
                $links{"$geB[0]!$geA[0]"}{"start_1!stop_2"}++;
            }else{
                $links{"$geA[0]!$geB[0]"}{$link}++;
            }
        }
        #print "$ke[0]\t$geA[0]\t$seqs{$geA[0]}{$ke[0]}[0][-2]\t$seqs{$geB[0]}{$ke[1]}[0][-1]\n";
        #print "Link is $link\n";
    }

    return;
}


################### Subroutines #####################

sub wait_for_file{
    my($fileName) = @_;

    # Using constants WAIT_TIME and NO_CHANGE_IN_FILE_FOR_X_SECONDS

    print "Wait_for_file subroutine invoked at ".localtime().", waiting for file $fileName\n";

    # Check for the finished file every WAIT_TIME minutes
    until (-e "$fileName") {
    print "$fileName does not exist yet at ".localtime()."\n";
    sleep WAIT_TIME;
    }

    print "File $fileName found at ".localtime()."\n";

    # make sure the file is no longer being modified (changing size)
    my$time = 0;
    my$size = -s $fileName;

    until ($time == NO_CHANGE_IN_FILE_FOR_X_SECONDS) {
        if ($size == -s $fileName) {
            $time++;
            sleep 1;
            if ($time%5 == 0) {
                print "No change in file size for $time seconds\n";
            }
        } else {
            $time = 0;
            $size = -s $fileName;
            print "The file size changed, sleeping for 1 minute\t";
            sleep 60;
            print"Waking up, try again\n";
        }
    }
    print "File $fileName exists and hasn't been modifed for at least 1 minute at time ".localtime()."\n\n";
    return $fileName;
}
