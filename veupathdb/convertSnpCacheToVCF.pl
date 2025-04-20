#!/usr/bin/env perl

use strict;
use warnings;
use lib './lib';

use ApiCommonData::Load::VariationFileReader;
# you will need:
# https://github.com/VEuPathDB/ApiCommonData/blob/master/Load/lib/perl/VariationFileReader.pm
# https://github.com/VEuPathDB/ApiCommonData/blob/master/Load/lib/perl/FileReader.pm
# https://github.com/VEuPathDB/ApiCommonData/blob/master/Load/lib/perl/SnpUtils.pm

use File::Temp qw/ tempfile tempdir /;
use Storable;

my @ordered_columns = ( # VCF output columns
  'CHROM',
  'POS',
  'ID',
  'REF',
  'ALT',
  'QUAL',
  'FILTER',
  'INFO',
  'FORMAT',
# Add column for each strain...
);


my $f = shift @ARGV;
my $ref_strain = shift @ARGV;

die "reference strain must be specified" unless $ref_strain;

my $vfr = ApiCommonData::Load::VariationFileReader->new($f, undef, qr/\t/);

my %all_strains;
my $strain_order = 1;

my @rows;

my $dir = tempdir(CLEANUP => 1); 

while( $vfr->hasNext() ){
## loop
  my ($fh, $filename) = tempfile('vcf_XXXXXX', DIR=> $dir);
  my $alleles = $vfr->nextSNP();
  
  ## Guess the reference allele
  ### !! really should get it from the reference genome fasta or db
  my $ref; # ref allele 
  foreach my $allele (@$alleles){
    if(($allele->{na_sequence_id} eq $allele->{ref_na_sequence_id}) &&
        ($allele->{matches_reference} &&
        $allele->{strain} eq $ref_strain)
      ){
      $ref = $allele;
      last;
    }
  }
  die "Cannot guess reference allele" unless $ref;
  my $rowdata = {
    CHROM => $ref->{sequence_source_id},
    POS => $ref->{location},
    ID => $ref->{snp_source_id},
    REF => $ref->{base},
    ALT => '',
    QUAL => '.',
    FILTER => '.',
    INFO => '.', 
    FORMAT => 'GT', 
  };
  my %alts; # accumulate, then sort descending for genotype
  foreach my $allele (@$alleles){
    ## count instances of base
    $alts{ $allele->{base} } ||= 0;
    $alts{ $allele->{base} }++;
    ## Validation check: if matches_reference, should be same base
    if ($allele->{matches_reference} ){
      if ($allele->{base} ne $ref->{base}){
        die "Invalid: matches_reference but base does not match";
      }
    }
  }
  # Normalize the allele mapping, 0 = ref, 1 = first alt, 2 = second alt, ...
  my $altnum = 0;
  my %genotype;
  foreach my $allele ( sort { $alts{$b} <=> $alts{$a} } keys %alts ){
    $genotype{ $allele } = $altnum;
    $altnum++;
  }
  # rowdata for each strain
  ## TODO: preserve the order of strains in output column order
  my %strain_alleles;
  foreach my $allele (@$alleles){
    if( defined($strain_alleles{ $allele->{strain} }) ){
      # only one per row!
      die "Duplicate strain";
    }
    $all_strains{ $allele->{strain} } ||= $strain_order++;
    $rowdata->{ $allele->{strain} } = $genotype{ $allele->{base} };
  }

  ## Interlude: Discard the reference base and write the ALT field
  delete %alts{ $ref->{base} };
  $rowdata->{ALT} = join(",", sort { $alts{$b} <=> $alts{$a} } keys %alts );
  $rowdata->{ALT} ||= '.';
  # push(@rows, $rowdata);
  store($rowdata, $filename);
  push(@rows, $filename);
}

unless(@rows){
  warn "No cache files written.";
  exit;
}

printf STDERR ("Writing output VCF...\n");
push( @ordered_columns, sort { $all_strains{$a} <=> $all_strains{$b} } keys %all_strains);

# header
printf( "%s\n", join("\n", 
  '##fileformat=VCFv4.2',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '#' . join("\t", @ordered_columns))
);
foreach my $filename (@rows){
  my $rowdata = retrieve($filename);
  printf( "%s\n", join("\t", map { defined($rowdata->{$_}) ? $rowdata->{$_} : '.' } @ordered_columns));
}  

__END__
