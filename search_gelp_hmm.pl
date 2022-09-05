#!/usr/bin/env perl

=pod

=head1 NAME

search_gelp_hmm.pl - search genome for GDSL genes and classify in orthogroups.

=head1 SYNOPSIS

    search_gelp_hmm.pl -d <hmm_dir> -f <fasta_file>

=head1 REQUIRES

Perl5
BioPerl
Bio::SearchIO::hmmer
forks

=head1 DESCRIPTION

Search genome(s) for GDSL genes and classify sequences in orthogroups.
Outputs (STDOUT) a 4 columns separated by tabulations:
Colunm 1: Matching HMM file name;
Colunm 2: Matching FASTA file name;
Colunm 3: Number of matching sequences;
Colunm 4: Coma separated list of matching sequence names.

=cut

use strict;
use warnings;
use forks;
use forks::shared;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SearchIO::hmmer;
use Getopt::Long;
use Pod::Usage;
use Carp qw (cluck confess croak);

++$|; #no buffering




# Script global constants
##########################

=pod

=head1 CONSTANTS

B<$EVALUE_INDEX>: (integer)

Index of best e-value in $gs_scan_results match sub-array.

B<$HMM_FILE_INDEX>: (integer)

Index of matched HMM file in $gs_scan_results match sub-array.

B<$HMM_NAME_INDEX>: (integer)

Index of matched HMM name in $gs_scan_results match sub-array.

B<$HMMSEARCH_COMMAND>: (string)

Command line to use to run hmmsearch. If hmmsearch is not in the PATH, provide
full path to the executable.

B<$DEFAULT_EVALUE>: (float)

Default e-value threshold.

=cut

our $EVALUE_INDEX = 0;
our $HMM_FILE_INDEX = 1;
our $HMM_NAME_INDEX = 2;
our $HMMSEARCH_COMMAND = 'hmmsearch';
# our $HMMSEARCH_COMMAND = 'module load hmmer/2.3.2; singularity exec /nfs/work/agap_id-bin/img/HMMER/2.3.2/hmmer.2.3.2.img hmmsearch';
our $DEFAULT_EVALUE = '1e-3';




# Script global variables
##########################

=pod

=head1 VARIABLES

B<$gs_evalue>: (float)

Selected HMM Search E-value threshold.

B<%gs_fasta_index>: (hash)

Contains numeric index (used in different arrays) of input FASTA files keyed by
their path.

B<%gs_hmm_index>: (hash)

Contains numeric index (used in different arrays) of HMM files keyed by
their path.

B<@gs_scan_stack>: (array)

Stack of HMM search (scans) to perform. Each value is an array of 2 values:
HMM file path to use as model and FASTA file path to scan.

B<@gs_scan_results>: (array)

First level of the array is indexed by fasta file names (through
%gs_fasta_index); Second level is a hash wich keys are sequence names; Third
level is an array containing best hit infos:
[$EVALUE_INDEX]: best e-value;
[$HMM_FILE_INDEX]: matching HMM file;
[$HMM_NAME_INDEX]: matching HMM name.
$gs_scan_results[<fasta file index>]->{<sequence name>} = [<evalue>, <HMM file>, <HMM name>]

B<@gs_hmm_best_matches>: (array)

First level is HMM file index and values are arrays of the name of best macthing
sequence and its e-value score for each HMM. This is used to detect sequences
best-matching more than one HMM profile.

=cut

my $gs_evalue : shared;
my %gs_fasta_index : shared;
my %gs_hmm_index : shared;
my @gs_scan_stack : shared;
my @gs_scan_results : shared;
my @gs_hmm_best_matches : shared;




# Script global functions
##########################

=pod

=head1 FUNCTIONS

=head2 run_scan_thread

B<Description>: Runs HMM search thread. The threads runs HMM search on the next
available set of HMM-FASTA file while there are some in the @gs_scan_stack
stack.

B<ArgsCount>: 0

B<Return>: nothing

=cut

sub run_scan_thread
{
    # Since this procedure is run from threads, modules should be reloaded.
    use strict;
    use warnings;
    use forks;
    use forks::shared;
    use Bio::SearchIO::hmmer;
    use Carp qw (cluck confess croak);

    my ($parameters) = @_;

    my $thread_id = threads->tid();
    # Foreach scan...
    while (@gs_scan_stack)
    {
        # Acquire a scan task using a lock system to prevent race conditions.
        my $current_scan;
        {
          lock(@gs_scan_stack);
          $current_scan = pop(@gs_scan_stack);
        }

        # Make sure we were able to obtain a task.
        if ($current_scan)
        {
            my ($hmm_file, $fasta_filename) = (@$current_scan);
            # Perform HMM search.
            my $hmm_result = `$HMMSEARCH_COMMAND $gs_evalue $hmm_file $fasta_filename`;
            # Check for errors.
            if (0 != $?)
            {
                warn "Failed to run HMM on $fasta_filename with $hmm_file profile: $!\n";
            }
            else
            {
                # Parse output result.
                my $fh;
                open($fh, '<', \$hmm_result);
                my $in = Bio::SearchIO->new(-format => 'hmmer', -fh => $fh);
                while (my $result = $in->next_result)
                {
                    while (my $hit = $result->next_hit)
                    {
                        # Acquire a lock on @gs_scan_results to update results.
                        {
                            lock(@gs_scan_results);
                            lock(@gs_hmm_best_matches);
                            # Process each hit.
                            while (my $hsp = $hit->next_hsp)
                            {
                                # Check if we already have a hit for that sequence.
                                if (exists($gs_scan_results[$gs_fasta_index{$fasta_filename}]->{$hit->name()}))
                                {
                                    # Yes we do, check if the new hit is better.
                                    if ($hsp->evalue() < $gs_scan_results[$gs_fasta_index{$fasta_filename}]->{$hit->name()}->[$EVALUE_INDEX])
                                    {
                                        # Better hit, replace previous one.
                                        $gs_scan_results[$gs_fasta_index{$fasta_filename}]->{$hit->name()} = &share([$hsp->evalue(), $hmm_file, $result->query_name(), ]);
                                    }
                                }
                                else
                                {
                                    # Nope, new hit.
                                    $gs_scan_results[$gs_fasta_index{$fasta_filename}]->{$hit->name()} = &share([$hsp->evalue(), $hmm_file, $result->query_name(), ]);
                                }
                                # Keep track of best matches.
                                if ((!defined($gs_hmm_best_matches[$gs_hmm_index{$hmm_file}]))
                                    || ($gs_hmm_best_matches[$gs_hmm_index{$hmm_file}][1] > $hsp->evalue()))
                                {
                                    $gs_hmm_best_matches[$gs_hmm_index{$hmm_file}] = &share([$hit->name(), $hsp->evalue()]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    # No more scans. We are done here.
    threads->exit();
}




# Script options
#################

=pod

=head1 OPTIONS

  check_genomes.pl
    [-help | -man]
    -f <FASTA_FILES>
    -d <HMM_DIRECTORY>
    -e <THRESOLD>
    -o <OUTPUT_DIRECTORY>
    -a
    -t <NB_THREADS>

=head2 Parameters

=over 4

=item B<-help>:

Prints a brief help message and exits.

=item B<-man>:

Prints the manual page and exits.

=item B<-f> or B<-fasta> or B<-genome> FASTA_FILES:

One or more amino acids (proteins) FASTA file(s) of genomes to search into.
Multiple FASTA files can be provided by using several "-f" parameters in the
command line.

=item B<-d> or B<-hmm-dir> HMM_DIRECTORY:

Path of the directory containg the HMM profile files to use for scanning.

=item B<-e> or B<-E> or B<-e-value> THRESOLD:

E-value threshold to use to discard poor HMM match. If not set, defaults to
1e-3.

=item B<-o> or B<-output> OUTPUT_DIRECTORY:

Optional. If set, OUTPUT_DIRECTORY must be a writable directory in wich FASTA
files of matching sequences will be generated. Each file will be named using the
matching HMM file name followd by a ".faa". If an output file already exists, it
will be overwritten unless the "-append" parameter is used. If you have more
than one input genome FASTA file, you should use "-append" otherwise, only the
matching sequences of the last FASTA will be output.

=item B<-a> or B<-append>:

Appends sequences to existing output FASTA files instead of overwritting the
files. Only useful when "-output" parameter is used.

=item B<-t> or B<-threads> NB_THREADS:

Specify a number of thread to use to run parallele scans. It is only usefull if
you have more than one HMM file and/or more than one FASTA file.
Default: 1 thread.

=back

=cut


# CODE START
#############

# options processing
my ($man, $help) = (0, 0);
my ($cmd, @fasta_filenames, $hmm_dir, $output_dir, $nb_slots, $append);
# parse options and print usage if there is a syntax error.
GetOptions('help|?'           => \$help,
           'man'              => \$man,
           'f|fasta|genome=s' => \@fasta_filenames,
           'd|hmm-dir=s'      => \$hmm_dir,
           'o|output:s'       => \$output_dir,
           'e|E|e-value=f'    => \$gs_evalue,
           't|threads:i'      => \$nb_slots,
           'a|append'         => \$append,
) or pod2usage(1);
if ($help) {pod2usage('-verbose' => 1, '-exitval' => 0);}
if ($man) {pod2usage('-verbose' => 2, '-exitval' => 0);}
if (!$hmm_dir || !@fasta_filenames) {pod2usage(1);}

# Check genome files.
my $i = 0;
foreach my $fasta_filename (@fasta_filenames)
{
    if (!-f $fasta_filename)
    {
        confess "ERROR: Genome FASTA file not found: $fasta_filename\n"
    }
    if (-z $fasta_filename)
    {
        confess "ERROR: Genome FASTA file is empty: $fasta_filename\n"
    }
    if (!-r $fasta_filename)
    {
        confess "ERROR: Genome FASTA file is not readable: $fasta_filename\n"
    }
    # Init %gs_fasta_index.
    $gs_fasta_index{$fasta_filename} = $i++;
    # Init @gs_scan_results first level.
    $gs_scan_results[$gs_fasta_index{$fasta_filename}] = &share({});
}

# Check if an output directory was specified.
if ($output_dir)
{
    if (!-d $output_dir)
    {
        confess "ERROR: Invalid output directory: $output_dir\n"
    }
    elsif (!-w $output_dir)
    {
        confess "ERROR: Output directory is not writable: $output_dir\n"
    }
    # Remove trailing slash.
    $output_dir =~ s#/+$##;
    if ('' eq $output_dir)
    {
        # Root directory was specified.
        $output_dir = '/.';
    }
}

# Check HMM directory.
if (!-d $hmm_dir)
{
    confess "ERROR: Invalid HMM directory: $hmm_dir\n"
}
# Trim trailing slashes.
$hmm_dir =~ s#/+$##;

# Get HMM files.
my @hmm_list = glob("$hmm_dir/*.hmm");
my @hmm_sequences;
# @result_table:
#   Level 1: fasta file index;
#   Level 2: HMM index;
#   Level 3: array of matching sequences.
my @result_table;
# Check HMM files.
$i = 0;
foreach my $hmm_file_path (@hmm_list)
{
    if (!-f $hmm_file_path || !-r $hmm_file_path)
    {
        confess "ERROR: Unable to access HMM file '$hmm_file_path'\n";
    }
    if (-z $hmm_file_path)
    {
        confess "ERROR: HMM file '$hmm_file_path' is empty\n";
    }
    $gs_hmm_index{$hmm_file_path} = $i++;
    $hmm_sequences[$gs_hmm_index{$hmm_file_path}] = [];
    
    # Prepare list of tasks.
    foreach my $fasta_filename (@fasta_filenames)
    {
        # Add scan task.
        push(@gs_scan_stack, [$hmm_file_path, $fasta_filename]);
        # Init @result_table first and second levels for later use.
        $result_table[$gs_fasta_index{$fasta_filename}] ||= [];
        $result_table[$gs_fasta_index{$fasta_filename}][$gs_hmm_index{$hmm_file_path}] = [];
    }
}

# Check e-value.
if (!$gs_evalue)
{
    # Make sure we can concatenate $gs_evalue.
    $gs_evalue = "-E $DEFAULT_EVALUE";
}
else
{
    # Check provided value.
    if (!($gs_evalue > 0) || ($gs_evalue !~ m/^[eE0-9\-+\.]$/)) {
        confess "ERROR: Invalid e-value '$gs_evalue'.\n";
    }
    # We add it as a parameter to hmmsearch command line.
    $gs_evalue = "-E $gs_evalue";
}

# Thread slots.
if (!$nb_slots || ($nb_slots < 1))
{
    # Use default.
    $nb_slots = 1;
}

my @threads_list = ();
# Launch threads.
warn "INFO: Searching using $nb_slots thread(s).\n";
for (my $slot_number = 0; $slot_number < $nb_slots; ++$slot_number)
{
    # Launch threads.
    push(@threads_list, threads->create(\&run_scan_thread));
}

# Wait for last threads to end.
foreach my $thread (@threads_list)
{
    $thread->join();
}

# Reorganize results by filling @result_table table.
foreach my $fasta_filename (@fasta_filenames)
{
    foreach my $sequence (keys(%{$gs_scan_results[$gs_fasta_index{$fasta_filename}]}))
    {
        my $hmm_file_path = $gs_scan_results[$gs_fasta_index{$fasta_filename}]->{$sequence}->[$HMM_FILE_INDEX];
        push(@{$result_table[$gs_fasta_index{$fasta_filename}][$gs_hmm_index{$hmm_file_path}]}, $sequence);
    }
}

# Check for multiple best matches.
my (%sequence_best_matches, %conflict_sequences);
foreach my $hmm_file_path (@hmm_list)
{
    my $hmm_file_name = $hmm_file_path;
    # Extract file name without path part.
    $hmm_file_name =~ s#.*/##;
  
    if (defined($gs_hmm_best_matches[$gs_hmm_index{$hmm_file_path}]))
    {
        my $sequence = $gs_hmm_best_matches[$gs_hmm_index{$hmm_file_path}][0];
        if (exists($sequence_best_matches{$sequence}))
        {
            $conflict_sequences{$sequence} ||= [$sequence_best_matches{$sequence}];
            push(@{$conflict_sequences{$sequence}}, $hmm_file_name);
        }
        else
        {
            $sequence_best_matches{$sequence} = $hmm_file_name;
        }
    }
}
# Reports conflicting sequences.
foreach my $conflict_sequence (keys(%conflict_sequences))
{
    # Note: conflicting sequences may be ok. Matched HMM listed here may also
    # report no matching sequence because the best sequence may have better
    # matched another HMM and will be reported there.
    warn "WARNING: sequence " . $conflict_sequence . " is the best match for: " . join(', ', @{$conflict_sequences{$conflict_sequence}}) . "\n";
}

# Output results.
foreach my $fasta_filepath (@fasta_filenames)
{
    my $fasta_filename = $fasta_filepath;
    # Extract file name without path part.
    $fasta_filename =~ s#.*/##;
    foreach my $hmm_file_path (@hmm_list)
    {
        my @sequence_list;
        if ($result_table[$gs_fasta_index{$fasta_filepath}][$gs_hmm_index{$hmm_file_path}])
        {
            @sequence_list = @{$result_table[$gs_fasta_index{$fasta_filepath}][$gs_hmm_index{$hmm_file_path}]};
        }
        # Extract file name without path part.
        my $hmm_file_name = $hmm_file_path;
        $hmm_file_name =~ s#.*/##;
        print
            $hmm_file_name . "\t"
            . $fasta_filename . "\t"
            . scalar(@sequence_list) . "\t"
            . (join(', ', @sequence_list) || '')
            . "\n";
        # Store sequences by HMM if we later need to output FASTA files.
        push(@{$hmm_sequences[$gs_hmm_index{$hmm_file_path}]}, @sequence_list);
    }
}

# Produce FASTA file by HMM profiles (optional).
if ($output_dir)
{
    # Load FASTAs to get original sequences.
    my %sequences;
    foreach my $fasta_filepath (@fasta_filenames)
    {
        my $fasta = Bio::SeqIO->new(
            -file   => $fasta_filepath,
            -format => 'fasta',
        );
        while (my $inseq = $fasta->next_seq)
        {
            if (exists($sequences{$inseq->display_id}))
            {
                confess "ERROR: Sequence '" . $inseq->display_id . "' found in more than one input FASTA file!\n";
            }
            $sequences{$inseq->display_id} = $inseq;
        }
    }

    # Output files.
    HMMFASTA: foreach my $hmm_file_path (@hmm_list)
    {
        my $hmm_file_name = $hmm_file_path;
        # Extract file name without path part.
        $hmm_file_name =~ s#.*/##;
        # Skip empty FASTA.
        if (!scalar(@{$hmm_sequences[$gs_hmm_index{$hmm_file_path}]}))
        {
            warn "INFO: no matches for $hmm_file_name; no FASTA generated.\n";
            next HMMFASTA;
        }
        # Generate a new file name.
        my $output_filename = $output_dir . '/' . $hmm_file_name . ".fa";
        if (!$append && -e $output_filename)
        {
            warn "WARNING: File '$output_filename' already exists and will be overwritten.\n";
        }
        my $fasta = Bio::SeqIO->new(
            -file   => ($append ? '>' : '') . ">$output_filename",
            -format => 'fasta',
        );
        foreach my $sequence (@{$hmm_sequences[$gs_hmm_index{$hmm_file_path}]})
        {
            if (exists($sequences{$sequence}))
            {
                $fasta->write_seq($sequences{$sequence});
            }
            else
            {
                warn "WARNING: Sequence '" . $sequence . "' has a different name in input FASTA files and could not be added to $output_filename.\n";
            }
        }
    }
}

exit(0);

__END__
# CODE END
###########

=pod

=head1 AUTHORS

Valentin GUIGNON (Bioversity), v.guignon@cgiar.org
Mathieu ROUARD (Bioversity), m.rouard@cgiar.org
Alberto CENCI (Bioversity), a.cenci@cgiar.org

=head1 VERSION

Version 1.0.0

Date 06/07/2022

=cut
