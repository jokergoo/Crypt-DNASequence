package Crypt::DNASequence;

use IO::Compress::Gzip qw(gzip);
use IO::Uncompress::Gunzip qw(gunzip);
use File::Temp qw(tempfile);
use strict;

our $VERSION = '0.2';

my $keys = [['00', '11', '01', '10'],
            ['00', '11', '10', '01'],
            ['00', '01', '11', '10'],
            ['00', '01', '10', '11'],
            ['00', '10', '11', '01'],
            ['00', '10', '01', '11'],
            ['11', '00', '01', '10'],
            ['11', '00', '10', '01'],
            ['11', '01', '00', '10'],
            ['11', '01', '10', '00'],
            ['11', '10', '00', '01'],
            ['11', '10', '01', '00'],
            ['01', '00', '10', '11'],
            ['01', '00', '11', '10'],
            ['01', '10', '00', '11'],
            ['01', '10', '11', '00'],
            ['01', '11', '00', '10'],
            ['01', '11', '10', '00'],
            ['10', '00', '01', '11'],
            ['10', '00', '11', '01'],
            ['10', '01', '00', '11'],
            ['10', '01', '11', '00'],
            ['10', '11', '00', '01'],
            ['10', '11', '01', '00']];

sub encrypt {
    my $class = shift;
    my $inputfile = shift;
	
	die "$inputfile is not a valid text file.\n" unless -T $inputfile;
	
	my ($tempfile_fh, $tempfile) = tempfile("XXXXXXXX", DIR => ".", SUFFIX => ".gz", UNLINK => 1);
	gzip($inputfile, $tempfile)
		or die "IO::Compress::Gzip failed: $IO::Compress::Gzip::GzipError\n";

    my $dict = initial_dict();
	
	my $gz_fh = $tempfile_fh;
	binmode($gz_fh);
	
	my $encryptedfile = "encrypted.fasta";
	
	open my $output_fh, ">", $encryptedfile or die $@;
	print $output_fh "  ";
	my $n_char = 2;
	
	my $first_letter;
	my $first_run_flag = 1;
	my $number_of_first_letter = 0;
	# read one byte each time
	while(read($gz_fh, my $buffer, 1)) {
		my $a = sprintf "%b", ord($buffer);
		my $binary_str = '0' x (8 - length($a)) . $a;
		
		my @nt = look_in_dict($binary_str, $dict);
		
		if($first_run_flag) {
			$first_letter = $nt[0];
			$first_run_flag = 0;
		}
		
		$number_of_first_letter += grep {$_ eq $first_letter} @nt;
		
		for (@nt) {
			$n_char %= 70;
			if($n_char == 0) {
				print $output_fh "\n";
			}
			
			print $output_fh $_;
			$n_char ++;
		}
	}

    my $key_index = $number_of_first_letter % scalar(@$keys);
    my @key_letters = map {$dict->{$_}} @{$keys->[$key_index]};
    
    for (@key_letters[2..3]) {
		$n_char %= 70;
		if($n_char == 0) {
			print $output_fh "\n";
		}
			
		print $output_fh $_;
		$n_char ++;
	}
	
	# move to the begining
	seek($output_fh, 0, 0);
	syswrite($output_fh, $key_letters[0], length($key_letters[0]));
	syswrite($output_fh, $key_letters[1], length($key_letters[1]));
    
    close($output_fh);
	
	print "Your encrypted sequence is in $encryptedfile\n";
}

sub decrypt {
    my $class = shift;
	my $inputfile = shift;
	
	die "$inputfile is not a valid file name.\n" unless -f $inputfile;
	
	my $inputfilesize = -s $inputfile;
	open my $input_fh, "<", $inputfile;
	my @key_letters;
	read($input_fh, $key_letters[0], 1);
	read($input_fh, $key_letters[1], 1);
	seek($input_fh, $inputfilesize - 2, 0);
	read($input_fh, $key_letters[2], 1);
	read($input_fh, $key_letters[3], 1);
	
    my $first_letter;
	seek($input_fh, 2, 0);
	read($input_fh, $first_letter, 1);
	
	my $number_of_first_letter = 1;
	for(my $i = 4; $i < $inputfilesize - 2; $i ++) {
		my $buffer;
		read($input_fh, $buffer, 1);
		if($buffer eq $first_letter) {
			$number_of_first_letter ++;
		}
	}
	
    my $key_index = $number_of_first_letter % scalar(@$keys);
    
    my $dict = {$key_letters[0] => $keys->[$key_index]->[0],
                $key_letters[1] => $keys->[$key_index]->[1],
                $key_letters[2] => $keys->[$key_index]->[2],
                $key_letters[3] => $keys->[$key_index]->[3]};
    
	my ($tempfile_fh, $tempfile) = tempfile("XXXXXXXX", DIR => ".", SUFFIX => ".gz", UNLINK => 1);
	
	my $gzfile = $tempfile;
	my $output_fh = $tempfile_fh;
	seek($input_fh, 2, 0);
	my $quater;
	my $i_nt = 0;
	for(my $i = 3; $i < $inputfilesize - 2; $i ++) {
		my $buffer;
		read($input_fh, $buffer, 1);
		
		next if($buffer eq "\n");
		
		$i_nt ++;
		$quater .= $dict->{$buffer};
		if($i_nt % 4 == 0) {
			print $output_fh chr(oct("0b$quater"));
			$quater = "";
		}
	}
	
	close($output_fh);
	
	my $decryptedfile = $gzfile;
	$decryptedfile = "decryptedfile";
	gunzip $gzfile, $decryptedfile;
	
	print "Your decrypted sequence is in $decryptedfile\n";
}

sub look_in_dict {
    my $b = shift;
    my $dict = shift;
    
    my @di = $b =~/(\d\d)/g;
    return (map {$dict->{$_}} @di);

}

sub initial_dict {
    
    if(rand() < 0.5) {
        return { _random_assign(['00', '11'], ['A', 'T']),
                 _random_assign(['10', '01'], ['C', 'G'])
                };
    }
    else {
        return { _random_assign(['00', '11'], ['C', 'G']),
                 _random_assign(['10', '01'], ['A', 'T'])
                };
    }
}

sub _random_assign {
    my $key = shift;
    my $value = shift;
    
    if(rand() < 0.5) {
        return ($key->[0] => $value->[0],
                $key->[1] => $value->[1]);
    }
    else {
        return ($key->[1] => $value->[0],
                $key->[0] => $value->[1]);
    }
}


__END__

=pod

=head1 NAME

Crypt::DNASequence - Encrypt and decrypt strings to DNA Sequences

=head1 SYNOPSIS

  use Crypt::DNASequence;
  
  Crypt::DNASequence->encrypt($input_filename);

  Crypt::DNASequence->decrypt($encrypted_file);

=head1 DESCRIPTION

The module is naiive and just for fun. It transforms text strings into DNA sequences. A DNA sequence
is composed of four nucleotides which are represented as A, T, C, G. If we transform 
"abcdefghijklmnopqistuvwxyzABCDEFGHIJKLMNOPQISTUVWXYZ", the corresponding sequence would be:

  TAGACCTGTCGGTGGGTGGTTCCCCGTATGAAGGGGGGGGCTATGAGTCTACAGACTGACAGGGGGAGTC
  AGCGAGTTAGCTAGTAAGCAAGTCCGCCCGTGCGCGCGTTCGCTCGTACGCACGTCGTCCGTTGCGCGGT
  CGGTCTGTTAGTCAGTTCTTCCTTTGTTCGACGTACAGACGTACATACGAACAAACGCCCACCCGGCCAG
  CCGTCCATCCGACCAACCGCGGCCGGTGCCAGGGCGGGCTGGTAGGCAGGTCTGCCTGTGTGCGGGGTGG
  GGCCCACGCCGCACAGCAGCAGGGGGGGGGGGGGGC

or

  CGAGTTCACTAACAAACAACCTTTTACGCAGGAAAAAAAATCGCAGACTCGTGAGTCAGTGAAAAAGACT
  GATAGACCGATCGACGGATGGACTTATTTACATATATACCTATCTACGTATGTACTACTTACCATATAAC
  TAACTCACCGACTGACCTCCTTCCCACCTAGTACGTGAGTACGTGCGTAGGTGGGTATTTGTTTAATTGA
  TTACTTGCTTAGTTGGTTATAATTAACATTGAAATAAATCAACGAATGAACTCATTCACACATAAAACAA
  AATTTGTATTATGTGATGATGAAAAAAAAAAAAAAT
  
The transformation is not unique due to a random mapping, but all the transformed sequences can be 
decrypted correctly to the origin string.

=head1 ALGORITHM

First, text file are compressed into gzip files when encrypting and DNA sequences
are first decrypted into gzip files and then uncompressed into normal text file.

The algorithm behind the module is simple. Two binary bits are used to represent a nucleotide such as '00' for A, '01' for C. 
If you have some knowledge of molecular biology, you would know that A only matches to T and C only matches to G.
So if '00' is choosen to be A, then '11' should be used to represent 'T'. In the module, the correspondence between binary bits
and nucleotides are applied randomly. The information of the correspondence dictionary is also stored in the finnal sequence.

Here is the procedure for encryption. 1. Split a string into a set of letters or charactors. 2. For each letter, convert to
its binary form and transform to ATCG every two bits using a randomly generated dictionary. The dictionary may looks like:

  $dict = { '00' => 'A',
            '11' => 'T',
            '01' => 'C',
            '10' => 'G' };

3. Join the A, T, G, C as a single sequence. 4. Find the first nucleotide of the sequence. 5, Find the number of the first nucleotide
in the sequence. 6. There is a database storing all arrangements of '00', '11', '01', '10'. 7. Calculate the index value from
the number of the first nucleotide by mod calculation. 8. Retrieve the arrangement with the index value, map them to the dictionary and get four nucleotides. E.g. the first nucleotide of the sequence is G. The number of G in the sequence is 40. The number of all arrangement
in the database is 24. Then we calculate the index value by 40 % 24 = 16. Then the 16th arrangement is retrieved and may looks like
['01', '11', '10', '00']. The four items in the array are mapped to the dictionary to be four nucleotides such as CTGA. Note this information
can be used in the decryption procedure.
9. Put the first two nucleotides at the begining of the sequence and the last two nucleotides at the end of the sequence. 10. That
is the finnal seuqence.

Here is the procedure for decryption. 1. Extract the first two and the last two nucleotides fromt the sequence. E.g. CT and GA. 
2. Count the number of the first nucleotide in the real sequence, e.g., 40 for G. 3. Use this number to calculate the index in the 
arrangement database, e.g., 16. 4. find the dictionary, i.e. a dictionary is generated from the 16th arrangement ['01', '11', '10', '00'] and
CTGA. 5. Translate the DNA sequence according the dictionary into binary bit form and finnaly to the orgin format.

=head2 Subroutines

=over 4

=item C<Crypt::DNASequence-E<gt>encrypt($input_file)>

encrypt the text file to DNA sequence

=item C<Crypt::DNASequence-E<gt>decrypt($encrypted_file)>

decrypt the DNA sequence to the origin text file.
               
=back

=head1 AUTHOR

Zuguang Gu E<lt>jokergoo@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2012 by Zuguang Gu

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.1 or,
at your option, any later version of Perl 5 you may have available.

=cut
