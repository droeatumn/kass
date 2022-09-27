#!/usr/bin/perl -w
# deep - a recursive directory traversal utility, by Neil Gunton
# This is free software, presented without any guarantee or warranty.

use strict 'vars';
use utf8;

# Get these modules from http://www.cpan.org if you don't have them already
use Getopt::Long;
use File::Glob ':glob';
use File::Path;

my $version = 'v2.11.3 2006-12-12';
my $copyright = 'Copyright 2002-2006 by Neil Gunton';
my $domain = 'http://www.neilgunton.com';

# Script for recursing over subdirectories and performing various operations on files
# which match a given pattern. Current commands are:
#  append    - recursively append a string to files matching given search pattern
#  chmod     - recursively chmod files matching a given pattern
#  clean     - recursively delete files ending in '~' or beginning & ending with '#' (usually left behind by editors)
#  delete    - recursively delete files matching a given pattern
#  do        - recursively do the given command on $filename, depending on optional perl condition
#  find      - recursively find a string in files matching a given pattern
#  output    - recursively outputs files matching filename pattern to stdout - useful for log analysis
#  prepend   - recursively prepend a string to files matching a given pattern
#  rename    - recursively rename files with a given extension to a different extension
#  replace   - recursively replace a string with another string in files matching a given pattern

# Just type the script name without any parameters (i.e. 'deep') to see syntax guide.

# The file pattern can be a list of patterns, separated by spaces, e.g. '*.c *.cpp'.
# If you do this, then you need to be sure to put quotes around this list on the command line.
# You should put quotes around any pattern on the command line which contains wildcards,
# otherwise the shell will expand this BEFORE passing it to Perl, which is not what we want,
# since the shell is only looking at the current directory, and we need the untouched wildcard for
# subdirectories.

## Syntax notes: 1. The vertical bar character is used to denote alternatives, e.g. '1|0' means '1 or 0'.
##               2. Command line options in [square brackets] are optional.

# Turn off screen output buffering, so we get stuff printed out immediately rather than in chunks
$| = 1;

# ------------------------------------------------------------------------------------------------------
# Config file routine
# Looks for .deep.cfg in . (current directory), then ~/ (user's home dir)
# Reads in file pattern string from the file, if found
# Also reads second line, list of directories to be excluded (separated by spaces)
# ------------------------------------------------------------------------------------------------------

sub get_config
{
    my ($pattern, $clean_pattern) = @_;

    $pattern = undef if ($pattern && $pattern =~ /^\-/);

    # Cater to both Linix and Windows environments
    my $home = $ENV{HOME} || $ENV{HOMEPATH};
    my $path = (-e './.deep.cfg') ? './.deep.cfg' : "$home/.deep.cfg";
    if (!(-e $path))
    {
	return ($pattern, undef, $clean_pattern);
    }
    else
    {
	open CONFIG, "< $path" or die "Could not open config file: $path: $!";

	# Main pattern is on first line, it can be overridden
	my $config_pattern = <CONFIG>;
	if ($config_pattern)
	{
	    chomp $config_pattern;
	}
	$pattern ||= $config_pattern;

	# Exclusion pattern is on optional 2nd line
	my $exclude = <CONFIG>;
	chomp $exclude;

	# Clean pattern is on optional 3rd line
	my $config_clean_pattern = <CONFIG>;
	if ($config_clean_pattern)
	{
	    chomp $config_clean_pattern;

	    # Only use it if it's not a blank line
	    if ($config_clean_pattern !~ /^\s*$/)
	    {
		$clean_pattern = $config_clean_pattern;
	    }
	}

	close CONFIG;

	return ($pattern, $exclude, $clean_pattern);
    }
}

# ------------------------------------------------------------------------------------------------------
# Main recursive directory traversal subroutines
# ------------------------------------------------------------------------------------------------------

sub build_hash
{
    my ($dir, $exclude_pattern) = @_;

    # Build list of files/directories to be excluded.
    # Since the exclusion list could include wildcards, we have to use the File::Glob routine and save the list for use later in the sub
    $exclude_pattern ||= '';
    my %exclude_hash;
    foreach my $exclusion_glob (split (/\s/, $exclude_pattern))
    {
	foreach my $file (File::Glob::bsd_glob ("$dir/$exclusion_glob"))
	{
	    # We want to keep note both of the file with directory on the front, and the filename by itself
	    $exclude_hash{$file} = 1;

	    # Backslash non word chars
	    my $dir_backslashed = $dir;
	    $dir_backslashed =~ s/([\W])/\\$1/g;

	    # Different OS's could use '/' (*nix) or '\' (Windows) to separate path components
	    $file =~ s/^${dir_backslashed}[\/\\]//i;
	    $exclude_hash{$file} = 1;
	}
    }

    return \%exclude_hash;
}

# ------------------------------------------------------------------------------------------------------

sub traverse
{
    my ($dir,            # The current directory which is to be processed
	$pattern,        # A pattern identifying files to be processed, e.g. '*.html *.epl'
	$process,        # What to process: 'f' = files, 'd' = dirs, 'fd' = both. Order and case not significant.
	$config_exclude, # A list of files/directories to be excluded from the search, e.g. '.svn .cvs' - from the config file, can be overridden by cmd line pattern
	$cmd_exclude,    # A list of files/directories to be excluded from the search, e.g. '.svn .cvs' - from command line, overrides command line pattern
	$subref,         # Subroutine reference, a callback function to process each file
	@sub_parms       # Parameter list to pass to the callback function
	) = @_;

    my $process_dirs = ($process =~ /d/i);
    my $process_files = ($process =~ /f/i);

    local *DIR;
    opendir (DIR, $dir) or die "Could not open directory $dir: $!";

    # Build hashrefs of the default config exclude pattern, and the command line exclusion pattern.
    # The config exclude pattern is overridden by the command line file pattern, which is in turn overridden by the command line exclusion pattern.
    my $config_exclude_hashref = $cmd_exclude ne '' ? {} : build_hash ($dir, $config_exclude);
    my $cmd_exclude_hashref = build_hash ($dir, $cmd_exclude);

    # Process files in this directory
    # Pattern consists of a potential list of patterns, separated by spaces.
    # First build list of patterns which user is explicitly requesting, so we can override the default config exclusion list later with the command line file pattern
    my @patterns = split (/\s/, $pattern);
    my %exclusion_override = ();
    foreach my $p (@patterns)
    {
	$exclusion_override{$p} = 1;
    }
 
    # Next we make a list of patterns, and then glob each of these
    foreach my $glob (@patterns)
    {
	# Iterate through the resulting list of files
	foreach my $file (File::Glob::bsd_glob ("$dir/$glob"))
	{
	    # Skip directories if $process_dirs is not set
	    my $exists  = -e $file;
	    my $is_dir  = -d $file;
	    my $is_file = ($exists && !$is_dir);
	    if ($exists &&
		(($is_dir && $process_dirs) || ($is_file && $process_files)) &&
		(!defined ($config_exclude_hashref->{$file}) || defined($exclusion_override{$glob})) &&
		!defined ($cmd_exclude_hashref->{$file}))
	    {
		# Call the handler subroutine
		$subref->($file, @sub_parms);
	    }
	}
    }

    # Now, recursively go down into subdirectories
    while (defined(my $file = readdir (DIR)))
    {
	# Only recurse on directories, which do not start with '.', and skip symbolic links
	if (-d "$dir/$file" &&
	    !(-l "$dir/$file") &&
	    ($file !~ /^\.{1,2}$/) &&
	    !defined ($config_exclude_hashref->{$file}) &&
	    !defined ($cmd_exclude_hashref->{$file}))
	{
	    traverse ("$dir/$file", $pattern, $process, $config_exclude, $cmd_exclude, $subref, @sub_parms);
	}
    }
}

# ------------------------------------------------------------------------------------------------------
# The command subroutines
# There are two subroutines for each command: cmd_<command> and callback_<command>
# cmd_command is the main subroutine which checks syntax and calls the main traversal routine
# callback_command is the subroutine which traverse() calls for each file
# ------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------
# COMMAND: append
# ------------------------------------------------------------------------------------------------------

sub cmd_append
{
    my ($string, $pattern) = @_;

    # Get the filename pattern
    my $config_exclude = '';
    ($pattern, $config_exclude) = get_config ($pattern);

    my $verbose = 1;
    my $cmd_exclude = '';
    GetOptions ('verbose=i' => \$verbose,
		'exclude=s' => \$cmd_exclude);

    if (!$pattern || !$string)
    {
	print qq{Syntax: deep append <string> [<file pattern>] [--verbose=1] [--exclude=<file pattern>]\n};
	print qq{  --verbose=1|0 - Whether files are printed as they are processed, default is 1\n};
	print qq{  --exclude=<file pattern> - Optional space separated file/directory patterns to be excluded (wildcards ok)\n};
	print qq{Append given string to files matching pattern\n};
	print qq{e.g. deep append 'some string' '*.txt'\n\n};
    }
    else
    {
	traverse ('.', $pattern, 'f', $config_exclude, $cmd_exclude, \&callback_append, ($string, $verbose));
    }
}

# -----------------------------------------------------------------------------------------------------

sub callback_append
{
    my ($filename, $string, $verbose) = @_;

    # Skip binary files
    if (-B $filename)
    {
	#print "Skipping binary file: $filename\n" if $verbose;
	return;
    }

    # Open up the input and temporary output files
    open (OUTFILE, ">> $filename") or die "Could not open $filename for appending: $!";

    # Write the string to the output file
    print OUTFILE $string;

    # Close file
    close (OUTFILE);
    print "Modified $filename\n";
}

# ------------------------------------------------------------------------------------------------------
# COMMAND: chmod
# ------------------------------------------------------------------------------------------------------

sub cmd_chmod
{
    my ($permissions, $pattern) = @_;

    my $config_exclude = '';
    my $dummy = '';
    ($dummy, $config_exclude) = get_config ($pattern);

    my $verbose = 1;
    my $process = 'fd';
    my $cmd_exclude = '';
    GetOptions ('verbose=i' => \$verbose,
		'process=s' => \$process,
		'exclude=s' => \$cmd_exclude);

    if (!$permissions || !$pattern)
    {
	print qq{Syntax: deep chmod <numeric permissions> <file pattern> [--verbose=1] [--process=fd] [--exclude=<file pattern>]\n};
	print qq{  --verbose=1|0 - Whether files are printed as they are processed, default is 1\n};
	print qq{  --process=f|d|fd - Process files, directories or both, default=fd (both)\n};
	print qq{  --exclude=<file pattern> - Optional space separated file/directory patterns to be excluded (wildcards ok)\n};
	print qq{Do chmod recursively on files/directories matching pattern\n};
	print qq{e.g. deep chmod 755 '*.html' --verbose=0\n};
    }
    else
    {
	traverse ('.', $pattern, $process, $config_exclude, $cmd_exclude, \&callback_chmod, ($permissions, $verbose));
    }
}

# ------------------------------------------------------------------------------------------------------

sub callback_chmod
{
    my ($filename, $permissions, $verbose) = @_;

    chmod (oct($permissions), $filename) or die "Could not chmod $filename: $!";
    print "chmod $permissions $filename\n" if $verbose;
}

# ------------------------------------------------------------------------------------------------------
# COMMAND: clean
# Deletes files ending in ~ (usually left behind by editors)
# ------------------------------------------------------------------------------------------------------

sub cmd_clean
{
    # Temp file pattern, can be overridden in config file
    my $dummy = '';
    my $config_exclude = '';
    my $clean_pattern = '*~ .*~ #*#';
    ($dummy, $config_exclude, $clean_pattern) = get_config (undef, $clean_pattern);

    my $verbose = 1;
    my $help = 0;
    my $cmd_exclude = '';
    GetOptions ('verbose=i' => \$verbose,
		'exclude=s' => \$cmd_exclude,
		'help' => \$help);

    if ($help)
    {
	print qq{Syntax: deep clean [--verbose=1|0] [--exclude=<file pattern>] [--help]\n};
	print qq{  --verbose=1|0 - Whether files are printed as they are processed, default is 1\n};
	print qq{  --exclude=<file pattern> - Optional space separated file/directory patterns to be excluded (wildcards ok)\n};
	print qq{  --help - print syntax, and the default pattern which will be used to delete files (override in config file)\n};
	print qq{Delete editor backup files (or whatever) which need to be cleaned up on a regular basis\n};
	print qq{e.g deep clean\n};
	print qq{Pattern: $clean_pattern\n};
    }
    else
    {
	# Recursively delete all files matching temp file pattern
	traverse ('.', $clean_pattern, 'f', $config_exclude, $cmd_exclude, \&callback_clean, ($verbose));
    }
}

# ------------------------------------------------------------------------------------------------------

sub callback_clean
{
    my ($filename, $verbose) = @_;

    unlink $filename or die "Could not delete $filename: $!";
    print "Deleted $filename\n" if $verbose;
}

# ------------------------------------------------------------------------------------------------------
# COMMAND: delete
# ------------------------------------------------------------------------------------------------------

sub cmd_delete
{
    my ($pattern) = @_;

    my $dummy = '';
    my $config_exclude = '';
    ($dummy, $config_exclude) = get_config ($pattern);

    # Get the options
    my $real = 0;
    my $verbose = 1;
    my $process = 'f';
    my $cmd_exclude = '';
    GetOptions ('real' => \$real,
		'verbose=i' => \$verbose,
		'process=s' => \$process,
		'exclude=s' => \$cmd_exclude);

    if (!$pattern)
    {
	print qq{Syntax: deep delete <file pattern> [--real] [--verbose=1|0] [--process=f|d|fd] [--exclude=<file pattern>]\n};
	print qq{  <pattern> - required file pattern (config file pattern is not used for this command)\n};
	print qq{  --real - Forces deletion 'for real'. Omit for dummy run, you'll be told what would be deleted.\n};
	print qq{  --verbose=1|0 - Whether files are printed as they are processed, default=1\n};
	print qq{  --process=f|d|fd - Process files, directories or both, default=f (files only)\n};
	print qq{  --exclude=<file pattern> - Optional space separated file/directory patterns to be excluded (wildcards ok)\n};
	print qq{Delete files and/or directories according to given pattern\n};
	print qq{e.g. deep delete '*.ico' --real\n};
    }
    else
    {
	traverse ('.', $pattern, $process, $config_exclude, $cmd_exclude, \&callback_delete, ($real, $verbose));
    }
}

# ------------------------------------------------------------------------------------------------------

sub callback_delete
{
    my ($filename, $real, $verbose) = @_;

    if (!$real)
    {
	if (-d $filename)
	{
	    print "!!! Would delete DIRECTORY: $filename\n" if $verbose;
	}
	else
	{
	    print "Would delete $filename\n" if $verbose;
	}
    }
    else
    {
	if (-d $filename)
	{
	    # Delete directory
	    rmtree ($filename) or die "Could not delete $filename: $!";
	}
	else
	{
	    # Delete file
	    unlink $filename or die "Could not delete $filename: $!";
	}
	print "Deleted $filename\n" if $verbose;
    }
}

# ------------------------------------------------------------------------------------------------------
# COMMAND: do
# ------------------------------------------------------------------------------------------------------

sub cmd_do
{
    my ($command, $pattern) = @_;

    my $config_exclude = '';
    ($pattern, $config_exclude) = get_config ($pattern);

    # Get options
    my $condition = '';
    my $verbose = 1;
    my $process = 'fd';
    my $cmd_exclude = '';
    GetOptions ('condition=s' => \$condition,
		'verbose=i' => \$verbose,
		'process=s' => \$process,
		'exclude=s' => \$cmd_exclude);

    if (!$command || !$pattern)
    {
	print qq{Syntax: deep do <command> [<file pattern>] [--condition=<perl condition>] [--verbose=1|0] [--process=f|d|fd] [--exclude=<file pattern>]\n};
	print qq{  <pattern> - optional file pattern, if not present then config file will be used\n};
	print qq{  --condition=<perl condition> - if present, then file is only processed if the condition evaluates to true.\n};
	print qq{                                 \$filename is set to path of current file being processed.\n};
	print qq{  --verbose=1|0 - Whether files are printed as they are processed, default=1\n};
	print qq{  --process=f|d|fd - Process files, directories or both, default=fd\n};
	print qq{  --exclude=<file pattern> - Optional space separated file/directory patterns to be excluded (wildcards ok)\n};
	print qq{Do given command recursively on files and directories\n};
	print qq{e.g. deep do 'gzip \$filename' '*.log' --condition='time - (stat(\$filename))[9] > 1000'\n};
    }
    else
    {
	traverse ('.', $pattern, $process, $config_exclude, $cmd_exclude, \&callback_do, ($command, $condition, $verbose));
    }
}

# ------------------------------------------------------------------------------------------------------

sub callback_do
{
    my ($filename, $command, $condition, $verbose) = @_;

    # Protect against '&', which makes shell think it's running something in the background
    $filename =~ s/\&/\\&/g;

    if ($condition eq '' || eval($condition))
    {
	# Make sure spaces in filenames are escaped
	$filename =~ s/(\s)/\\$1/g;

	# Replace variable with filename in command
	$command =~ s/\$filename/$filename/g;

	# Do command
	print "$command\n" if $verbose;
	eval ("system \'$command\'") == 0 or die "$command failed: $?";
    }
}

# ------------------------------------------------------------------------------------------------------
# COMMAND: find
# ------------------------------------------------------------------------------------------------------

sub cmd_find
{
    my ($search_string, $pattern) = @_;

    my $config_exclude = '';
    ($pattern, $config_exclude) = get_config ($pattern);

    # Get the options
    my $case = 1;
    my $literal = 1;
    my $delete = undef;
    my $verbose = 1;
    my $cmd_exclude = '';
    GetOptions ('case=i'    => \$case,
		'literal=i' => \$literal,
		'delete'    => \$delete,
		'verbose=i' => \$verbose,
		'exclude=s' => \$cmd_exclude);

    if (!$pattern || !$search_string)
    {
	print qq{Syntax: deep find <search string> [<file pattern>] [--case=1|0] [--delete] [--verbose=1|0] [--exclude=<file pattern>]\n};
	print qq{  --case=1      case sensitive search (default)\n};
	print qq{  --case=0      case insensitive search\n};
	print qq{  --delete      delete files which contain a match\n};
	print qq{  --literal=1   backslash special pattern-match characters, i.e. treat string literally (default)\n};
	print qq{  --literal=0   do not backslash special pattern-match characters, i.e. allow regular expressions\n};
	print qq{  --exclude=<file pattern> - Optional space separated file/directory patterns to be excluded (wildcards ok)\n};
	print qq{Search for string in text files\n};
	print qq{e.g. deep find 'hooty' '*.html *.epl' --case=0 --delete\n};
    }
    else
    {
	if ($literal)
	{
	    # Backslash non word chars
	    $search_string =~ s/([\W])/\\$1/g;
	}

	traverse ('.', $pattern, 'f', $config_exclude, $cmd_exclude, \&callback_find, ($search_string, $case, $literal, $delete, $verbose));
    }
}

# ------------------------------------------------------------------------------------------------------

sub callback_find
{
    my ($filename, $search_string, $case, $literal, $delete, $verbose) = @_;

    # Skip binary files
    if (-B $filename)
    {
	#print "Skipping binary file: $filename\n" if $verbose;
	return;
    }

    # Open up the input file
#droe    open (INFILE, "< $filename") or die "Could not open $filename for reading: $!";
    open (INFILE, "< $filename") or print "Could not open $filename for reading: $!";
    
    # Keep note of whether we have found the string in this file yet
    my $found = 0;

    # Read in all the lines of the input file
    while (my $line = <INFILE>)
    {
	# Search for the string in the current line
	my $result = $case ? ($line =~ /$search_string/) : ($line =~ /$search_string/i);

	# Print out the filename, if this is the first occurance
	if ($result && !$found)
	{
	    print "$filename\n";
	}
	
	# Take note of whether we found anything yet - so we can tell the first occurrence
	$found = $found ? $found : $result;

	# Print out the string if it was found
	if ($result)
	{
	    # Get rid of any leading tabs or spaces
	    $line =~ s/^\s*//;
	    print "\t$.\t$line";
	}
    }

    # If this is finddel or finddeli, then delete the file
    if ($found && $delete)
    {
	unlink $filename or die "Could not delete $filename: $!";
	print "Deleted $filename\n" if $verbose;
    }

    close (INFILE);
}

# ------------------------------------------------------------------------------------------------------
# COMMAND: output
# ------------------------------------------------------------------------------------------------------

sub cmd_output
{
    my ($pattern) = @_;

    # Get the filename pattern
    my $config_exclude = '';
    ($pattern, $config_exclude) = get_config ($pattern);

    my $verbose = 1;
    my $cmd_exclude = '';
    GetOptions ('verbose=i' => \$verbose,
		'exclude=s' => \$cmd_exclude);

    if (!$pattern)
    {
	print qq{Syntax: deep output [<file pattern>] [--verbose=1|0] [--exclude=<file pattern>]\n};
	print qq{  <pattern> - optional file pattern, if not present then config file will be used\n};
	print qq{  --verbose=1|0 - Whether files are printed as they are processed, default=1\n};
	print qq{  --exclude=<file pattern> - Optional space separated file/directory patterns to be excluded (wildcards ok)\n};
	print qq{Output to stdout contents of matching text files\n};
	print qq{e.g. deep output '*.txt'\n};
    }
    else
    {
	traverse ('.', $pattern, 'f', $config_exclude, $cmd_exclude, \&callback_output, ($verbose));
    }
}

# -----------------------------------------------------------------------------------------------------

sub callback_output
{
    my ($filename, $verbose) = @_;

    # Skip binary files
    if (-B $filename)
    {
	#print "Skipping binary file: $filename\n" if $verbose;
	return;
    }

    # Open up the input file
    open (INFILE, "< $filename") or die "Could not open $filename for reading: $!";

    # Read in all the lines of the input file and write them directly to output
    while (my $line = <INFILE>)
    {
	# Output the line to the temp file
	print $line;
    }
    close (INFILE);
}

# ------------------------------------------------------------------------------------------------------
# COMMAND: prepend
# ------------------------------------------------------------------------------------------------------

sub cmd_prepend
{
    my ($string, $pattern) = @_;

    # Get the filename pattern
    my $config_exclude = '';
    ($pattern, $config_exclude) = get_config ($pattern);

    my $verbose = 1;
    my $cmd_exclude = '';
    GetOptions ('verbose=i' => \$verbose,
		'exclude=s' => \$cmd_exclude);

    if (!$pattern || !$string)
    {
	print qq{Syntax: deep prepend <string> [<file pattern>] [--verbose=1|0] [--exclude=<file pattern>]\n};
	print qq{  --verbose=1|0 - Whether files are printed as they are processed, default=1\n};
	print qq{  --exclude=<file pattern> - Optional space separated file/directory patterns to be excluded (wildcards ok)\n};
	print qq{Prepend text to start of matching text files\n};
	print qq{e.g. deep prepend 'some string' '*.txt'\n};
    }
    else
    {
	traverse ('.', $pattern, 'f', $config_exclude, $cmd_exclude, \&callback_prepend, ($string, $verbose));
    }
}

# -----------------------------------------------------------------------------------------------------

sub callback_prepend
{
    my ($filename, $string, $verbose) = @_;

    # Skip binary files
    if (-B $filename)
    {
	#print "Skipping binary file: $filename\n" if $verbose;
	return;
    }

    # Open up the input and temporary output files
    open (INFILE, "< $filename") or die "Could not open $filename for reading: $!";
    open (OUTFILE, "> $filename.tmp") or die "Could not open $filename.tmp for writing: $!";

    # Write the string to the output file
    print OUTFILE $string;

    # Read in all the lines of the input file and write them directly to output
    while (my $line = <INFILE>)
    {
	# Output the line to the temp file
	print OUTFILE $line;
    }
    close (INFILE);
    close (OUTFILE);

    # Get the uid and gid of the old file
    my ($mode, $uid, $gid) = (stat ($filename))[2,4,5];

    # Set the ownership of the new file to be same as old file
    chown ($uid, $gid, "$filename.tmp") == 1 or die "chown failed for $filename.tmp";
    chmod ($mode, "$filename.tmp") == 1 or die "chmod failed for $filename.tmp";

    # Rename the temp file to replace the old file
    rename ("$filename.tmp", $filename) or die "Could not rename $filename.tmp to $filename: $!";
    print "Modified $filename\n" if $verbose;
}

# ------------------------------------------------------------------------------------------------------
# COMMAND: rename
# ------------------------------------------------------------------------------------------------------

sub cmd_rename
{
    my ($old_ext, $new_ext) = @_;

    my $pattern = '';
    my $config_exclude = '';
    ($pattern, $config_exclude) = get_config();

    my $verbose = 1;
    my $process = 'fd';
    my $svn = 0;
    my $cmd_exclude = '';
    GetOptions ('verbose=i' => \$verbose,
		'process=s' => \$process,
		'exclude=s' => \$cmd_exclude,
		'svn=i'     => \$svn);

    if (!$old_ext || !$new_ext)
    {
	print qq{Syntax: deep rename <file extension> <new extension> [--verbose=1|0] [--process=f|d|fd] [--exclude=<file pattern>] [--svn=1|0]\n};
	print qq{  --verbose=1|0 - Whether files are printed as they are processed, default=1\n};
	print qq{  --process=f|d|fd - Process files, directories or both, default=fd (both)\n};
	print qq{  --exclude=<file pattern> - Optional space separated file/directory patterns to be excluded (wildcards ok)\n};
	print qq{  --svn=[1|0] - Optional - do svn rename rather than normal rename - for svn directories\n};
	print qq{Rename filenames ending in matching pattern\n};
	print qq{e.g. deep rename JPEG jpg\n};
    }
    else
    {
	traverse ('.', "*$old_ext", $process, $config_exclude, $cmd_exclude, \&callback_rename, ($old_ext, $new_ext, $verbose, $svn));
    }
}

# ------------------------------------------------------------------------------------------------------

sub callback_rename
{
    my ($filename, $old_ext, $new_ext, $verbose, $svn) = @_;

    if ($filename =~ /$old_ext$/)
    {
	my $new_filename = $filename;
	$new_filename =~ s/$old_ext$/$new_ext/;
	if ($svn)
	{
	    system ("svn rename $filename $new_filename");
	}
	else
	{
	    rename ($filename, $new_filename) or die "Could not rename $filename to $new_filename: $!";
	}
	print "Renamed $filename => $new_filename\n" if $verbose;
    }
}

# ------------------------------------------------------------------------------------------------------
# COMMAND: replace
# ------------------------------------------------------------------------------------------------------

sub cmd_replace
{
    my ($from, $to, $pattern) = @_;

    my $config_exclude = '';
    ($pattern, $config_exclude) = get_config ($pattern);

    # Get the options
    my $case = 1;
    my $literal = 1;
    my $raw = 1;
    my $verbose = 1;
    my $cmd_exclude = '';
    GetOptions ('case=i'  => \$case,
		'literal=i' => \$literal,
		'raw=i' => \$raw,
		'verbose=i' => \$verbose,
		'exclude=s' => \$cmd_exclude);

    if (!$pattern || !$from)
    {
	print qq{Syntax: deep replace <search string> <replace string> [<file pattern>] [--case=1|0] [--literal=1|0] [--raw=1|0] [--verbose=1|0] [--exclude=<file pattern>]\n};
	print qq{  --case=1      case sensitive search (default)\n};
	print qq{  --case=0      case insensitive search\n};
	print qq{  --literal=1   backslash special pattern-match characters, i.e. treat string literally (default)\n};
	print qq{  --literal=0   do not backslash special pattern-match characters, i.e. allow regular expressions\n};
	print qq{  --raw=1       treat input and output streams as "raw" data (default)\n};
	print qq{  --raw=0       do not treat input and output streams as "raw" data (for when UTF-8 support works reliably)\n};
	print qq{  --exclude=<file pattern> - Optional space separated file/directory patterns to be excluded (wildcards ok)\n};
	print qq{Replace text in matching files\n};
	print qq{e.g. deep replace 'hooty' 'blowfish' '*.html' --case=0\n};
    }
    else
    {
	if ($literal)
	{
	    # Backslash non word chars
	    $from =~ s/([\W])/\\$1/g;
	}

	traverse ('.', $pattern, 'f', $config_exclude, $cmd_exclude, \&callback_replace, ($raw, $case, $literal, $from, $to, $verbose));
    }
}

# -----------------------------------------------------------------------------------------------------

sub callback_replace
{
    my ($filename, $raw, $case, $literal, $from, $to, $verbose) = @_;

    # Skip binary files
    if (-B $filename)
    {
	# Get the options
	#print "Skipping binary file: $filename\n" if $verbose;
	return;
    }

    # Open up the input and temporary output files
    open (INFILE, "< $filename") or die "Could not open $filename for reading: $!";
    binmode (INFILE, ":raw") if $raw;
    open (OUTFILE, "> $filename.tmp") or die "Could not open $filename.tmp for writing: $!";
    binmode (OUTFILE, ":raw") if $raw;
    
    # Keep note of whether we have modified the file
    my $changed = 0;

    # Read in all the lines of the input file
    while (my $line = <INFILE>)
    {
	# Do the replace operation on the current line
	my $result = $case ? ($line =~ s/$from/$to/g) : ($line =~ s/$from/$to/gi);
	
	# Take note of whether the file has been changed
	$changed = $changed ? $changed : $result;

	# Output the (possibly) modified line to the temp file
	print OUTFILE $line;
    }
    close (INFILE);
    close (OUTFILE);

    if ($changed)
    {
	# Get the uid and gid of the old file
	my ($mode, $uid, $gid) = (stat ($filename))[2,4,5];

	# Set the ownership and permissions of the new file to be same as old file
	chown ($uid, $gid, "$filename.tmp") == 1 or die "chown failed for $filename.tmp";
	chmod ($mode, "$filename.tmp") == 1 or die "chmod failed for $filename.tmp";

	# Delete the old file, and rename the new one
	rename ("$filename.tmp", $filename) or die "Could not rename $filename.tmp to $filename: $!";
	print "Modified $filename\n" if $verbose;
    }
    else
    {
	# Remove the temp file
	unlink "$filename.tmp" or die "Could not remove $filename.tmp: $!";
    }
}

# -----------------------------------------------------------------------------------------------------

sub cmd_help
{
do_syntax();
print <<ENDHELP

FILE PATTERNS:

Can be a list of wildcard patterns, separated by spaces, e.g. '*.c
*.cpp'.  Always put single quotes around wildcard file patterns,
e.g. '*.html' rather than just *.html.  This is because some shells
automatically expand patterns before passing them to Perl.

CONFIG FILE:

The script first looks for .deep.cfg in the current directory and then
in the user's home directory. If found then this file can contain up
to three lines of plain text. The first line is the default filename
pattern for the commands which can use this ('append', 'find',
'output', 'prepend' and 'replace'). This allows you to omit the
filename pattern on the command line for these commands. An optional
second line may contain a space-separated list of file/directory
patterns (wildcards are ok) that should NOT be traversed or
processed. The third line, if present, allows you to override the
pattern list for the 'clean' command. For example:

*.html *.epl *.pm *.h *.hpp *.cpp *.c
.svn
*~ .*~ #*# ~\$*

If you want to specify line 3 but not line 2, then just leave line 2
empty.

ENDHELP
}

# -----------------------------------------------------------------------------------------------------

sub do_syntax
{
print <<ENDSYNTAX

deep utility $version
$copyright: $domain

Syntax:
    deep append <string> [<file pattern>] [--verbose=1|0] [--exclude=<file pattern>]
    deep chmod <numeric permissions> <file pattern> [--verbose=1|0] [--process=f|d|fd] [--exclude=<file pattern>]
    deep clean [--verbose=1|0] [--exclude=<file pattern>] [--help]
    deep delete <file pattern> [--real] [--verbose=1|0] [--process=f|d|fd] [--exclude=<file pattern>]
    deep do <command> [<file pattern>] [--condition=<perl condition>] [--verbose=1|0] [--process=f|d|fd] [--exclude=<file pattern>]
    deep find <search string> [<file pattern>] [--case=1|0] [--literal=1|0] [--delete] [--verbose=1|0] [--exclude=<file pattern>]
    deep help
    deep output [<file pattern>] [--verbose=1|0] [--exclude=<file pattern>]
    deep prepend <string> [<file pattern>] [--verbose=1|0] [--exclude=<file pattern>]
    deep rename <file extension> <new extension> [--verbose=1|0] [--process=f|d|fd] [--exclude=<file pattern>] [--svn=1|0]
    deep replace <search string> <replace string> [<file pattern>] [--case=1|0] [--literal=1|0] [--raw=1|0] [--verbose=1|0] [--exclude=<file pattern>]

Type a command by itself for detailed syntax on that command, e.g. 'deep replace'.

'deep help' for more help

ENDSYNTAX
}

# =====================================================================================================
# Main code
# =====================================================================================================

{
    # This is a hash which maps the command name to its subroutine
    my %commands = (
		    'append'   => \&cmd_append,
		    'chmod'    => \&cmd_chmod,
		    'clean'    => \&cmd_clean,
		    'delete'   => \&cmd_delete,
		    'do'       => \&cmd_do,
		    'find'     => \&cmd_find,
		    'help'     => \&cmd_help,
		    'output'   => \&cmd_output,
		    'prepend'  => \&cmd_prepend,
		    'rename'   => \&cmd_rename,
		    'replace'  => \&cmd_replace,
		    );

    # There should always be a command
    my $command = shift;

    # Do the command, or print an error message
    if (!$command || !$commands{$command})
    {
	do_syntax();
    }
    else
    {
	$commands{$command}->(@ARGV);
    }
}

1;

# ------------------------------------------------------------------------------------------------------
# End of file
# ------------------------------------------------------------------------------------------------------
