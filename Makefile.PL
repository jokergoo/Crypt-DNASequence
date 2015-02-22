use ExtUtils::MakeMaker;

WriteMakefile(
	NAME          => 'Crypt::DNASequence',
	VERSION_FROM  => 'lib/Crypt/DNASequence.pm',
	PREREQ_PM     => {
		'Test::More'   => 0,
		'IO::Compress::Gzip' => 0,
		'IO::Uncompress::Gunzip' => 0,
		'File::Temp' => 0,
	},
	AUTHOR        => 'Zuguang Gu <jokergoo@gmail.com>',
	ABSTRACT      => 'Encrypt and decrypt string to DNA sequence',
	LICENSE       => 'perl',
);
