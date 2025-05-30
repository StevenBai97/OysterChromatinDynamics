use strict;
use warnings;

# Variables 
my $organism = 'Cni';
my $pattern = 'Cni_';
my $nodes = 'Cni_nodes.txt';
my $outfile = 'Cni_paralogs.txt';
my $countsfile = 'Cni_paralog_counts.txt';
my $no_paralogs = 'Cni_missed.txt';

# Open outfiles
open OUT, ">$outfile" or die $!;
print OUT "gene1	gene2	duplication_node	Orthogroup	Species Tree Node	Gene Tree Node	Support	Type\n";

open OUTNO, ">$no_paralogs" or die $!;
print OUTNO "gene1	gene2	duplication_node	Orthogroup	Species Tree Node	Gene Tree Node	Support	Type	Reason for exclusion\n";

open CT, ">$countsfile" or die $!;
print CT "duplication_node	name	genes\n";

# Read the nodes file and make a hash
open FH, $nodes or die $!;
my @nodes = <FH>;
shift @nodes;
my @node_order;
my %Cni_nodes;

foreach (@nodes) {
    my @line = split "\t", $_;
    map { $_ =~ s/\n//g } @line;
    
    $Cni_nodes{$line[0]} = $line[1];
    push @node_order, $line[0];
}

# Read paralog file
open FH2, "../Gene_Duplication_Events/Duplications.tsv" or die $!;
my %uniqueNodes;

foreach (<FH2>) {
    my @line = split "\t", $_;
    map { $_ =~ s/\n//g } @line;
    
    if (exists $Cni_nodes{$line[1]}) {
        
        if (($line[5] =~ /$organism/g) && ($line[6] =~ /$organism/g)) {
            my @gene1 = split ',', $line[5];
            my @gene2 = split ',', $line[6];
            map { $_ =~ s/\n|\s+//g } @gene1;
            map { $_ =~ s/\n|\s+//g } @gene2;

            my @Cni1 = grep { $_ =~ /$organism/ } @gene1;
            map { $_ =~ s/$pattern//g } @Cni1;

            my @Cni2 = grep { $_ =~ /$organism/ } @gene2;
            map { $_ =~ s/$pattern//g } @Cni2;

            foreach my $g1 (@Cni1) {
                foreach my $g2 (@Cni2) {
                    print OUT "$g1	$g2	$Cni_nodes{$line[1]}	$line[0]	$line[1]	$line[2]	$line[3]	$line[4]\n";
                    $uniqueNodes{$line[1]} += 1;
                }
            }		
        }
        elsif (($line[5] =~ /$organism/g) && ($line[6] !~ /$organism/g)) {
            print OUTNO join("\t", @line[0..4]), "\tgene2 does not have Cni gene\n";
        }
        elsif (($line[5] !~ /$organism/g) && ($line[6] =~ /$organism/g)) {
            print OUTNO join("\t", @line[0..4]), "\tgene1 does not have Cni gene\n";
        }
        else {
            print OUTNO join("\t", @line[0..4]), "\tboth genes do not have Cni gene\n";
        }
    }
}

# Print node counts in order
my %checked;
foreach (@node_order) {
    print CT "$_\t$Cni_nodes{$_}\t$uniqueNodes{$_}\n";
    if (not exists $uniqueNodes{$_}) { die "$_ not found in nodes\n" }
    $checked{$_} = '';
}

foreach (keys %uniqueNodes) {
    if (not exists $checked{$_}) {
        print CT "\t$_\t$Cni_nodes{$_}\t$uniqueNodes{$_}\n";
    }
}

print `date`;
