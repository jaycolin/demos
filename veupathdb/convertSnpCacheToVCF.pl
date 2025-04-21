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
=head1 NAME

convertSnpCacheToVCF.pl - Convert SNP cache data into VCF (Variant Call Format) file

=head1 SYNOPSIS

  convertSnpCacheToVCF.pl <input_file> <reference_strain>
  
  Example:
  convertSnpCacheToVCF.pl snp_data.txt reference_strain_name > output.vcf

=head1 DESCRIPTION

This script reads SNP (Single Nucleotide Polymorphism) data from a cache file and converts it into a VCF (Variant Call Format) file. The VCF file is a widely used format for storing gene sequence variations.

=head1 REQUIREMENTS

The script requires the following dependencies, which can be found in the VEuPathDB GitHub repository:

=over 4

=item * L<ApiCommonData::Load::VariationFileReader> - Reads variation files.

=item * L<ApiCommonData::Load::FileReader> - Base file reader library.

=item * L<ApiCommonData::Load::SnpUtils> - Utilities for SNP data processing.

=back

=head1 ARGUMENTS

=over 4

=item B<input_file>

The input SNP data file.

=item B<reference_strain>

The name of the reference strain, which must be specified.

=back

=head1 OUTPUT

The script outputs a VCF file to STDOUT. The VCF file contains the following columns:

=over 4

=item * CHROM - The chromosome or scaffold name.

=item * POS - The position of the SNP.

=item * ID - The SNP identifier.

=item * REF - The reference allele.

=item * ALT - The alternative allele(s).

=item * QUAL - Placeholder for quality score (set to ".").

=item * FILTER - Placeholder for filter status (set to ".").

=item * INFO - Placeholder for additional information (set to ".").

=item * FORMAT - The genotype format (set to "GT").

=item * Additional columns for strains, with genotype data.

=back

=head1 DETAILS

The script processes the input SNP data in the following steps:

=over 4

=item 1. Reads SNP data using L<ApiCommonData::Load::VariationFileReader>.

=item 2. Identifies the reference allele for each SNP. If the reference allele cannot be determined, the script terminates with an error.

=item 3. Generates genotype data for each strain, normalizing allele mappings (0=reference, 1=first alternative, etc.).

=item 4. Outputs the SNP data in VCF format with the required columns and headers.

=back

=head1 TEMPORARY FILES

The script uses temporary files to store intermediate SNP data and removes them upon completion.

=head1 LIMITATIONS

=over 4

=item * The script attempts to guess the reference allele. Ideally, the reference allele should be obtained from a reference genome FASTA or database.

=item * The script assumes unique strain entries per SNP. Duplicate strains in the input data will result in an error.

=item * The order of strains in the output VCF is determined by the order of their appearance in the input data.

=back

=head1 EXAMPLES

  # Convert SNP data to VCF format
  perl convertSnpCacheToVCF.pl my_snp_data.txt my_reference_strain > output.vcf

=head1 AUTHOR

This script is part of the VEuPathDB project.

=head1 LICENSE

This script is distributed under the terms of the VEuPathDB license.

=cut
